# derived from minimap2.py in http://github.com/PoonLab/covizu

import subprocess
import argparse
import re
import sys
import json
from datetime import date, timedelta
from scipy.stats import poisson
from epiweeks import Week


def convert_fasta(handle):
    """
    Parse FASTA file as a list of header, sequence list objects
    :param handle:  open file stream
    :return:  List of [header, sequence] records
    """
    result = []
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    result.append([h, sequence])  # handle last entry
    return result


def load_filter(vcf_file):
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.

    :param features:  list, return object from import_json()
    :param vcf_file:  str, path to VCF file
    :return:  list, filtered features
    """
    vcf = open(vcf_file)
    mask = {}
    for line in vcf.readlines():
        if line.startswith('#'):
            continue
        _, pos, _, ref, alt, _, filt, info = line.strip().split()
        if filt == 'mask':
            mask.update({int(pos)-1: {  # convert to 0-index
                'ref': ref, 'alt': alt, 'info': info}
            })
    return mask


def apply_cigar(seq, rpos, cigar):
    """
    Use CIGAR to pad sequence with gaps as required to
    align to reference.  Adapted from http://github.com/cfe-lab/MiCall
    """
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
    if not is_valid:
        raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
    aligned = '-'*rpos
    left = 0
    for length, operation in tokens:
        length = int(length)
        if operation in 'M=X':
            aligned += seq[left:(left+length)]
            left += length
        elif operation == 'D':
            aligned += '-'*length
        elif operation in 'SI':
            left += length  # soft clip

    return aligned


def minimap2(fasta, ref, path='minimap2', nthread=3, minlen=29000):
    """
    Wrapper function for minimap2.

    :param fasta:  str, path to FASTA with query sequences
    :param ref:  str, path to FASTA with reference sequence(s)
    :param path:  str, path to binary executable
    :param nthread:  int, number of threads for parallel execution of minimap2
    :param minlen:  int, filter genomes below minimum length; to accept all, set to 0.

    :yield:  query sequence name, reference index, CIGAR and original
             sequence
    """
    p = subprocess.Popen([path, '-t', str(nthread), '-a', '--eqx', ref, fasta],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in map(lambda x: x.decode('utf-8'), p.stdout):
        if line.startswith('@'):
            continue
        qname, flag, rname, rpos, _, cigar, _, _, _, seq = \
            line.strip().split()[:10]
        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            continue

        if len(seq) < minlen:
            # reject sequence that is too short
            continue

        # validate CIGAR string
        is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
        if not is_valid:
            raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

        rpos = int(rpos) - 1  # convert to 0-index
        yield qname, rpos, cigar, seq


def output_fasta(iter, outfile, reflen=0):
    """
    Stream output from minimap2 into FASTA file
    of aligned sequences.  CIGAR parsing code adapted from
    http://github.com/cfe-lab/MiCall

    :param iter:  generator from minimap2()
    :param outfile:  open file stream in write mode
    :param reflen:  int, length of reference genome to pad sequences;
                    defaults to no padding.
    """
    for qname, rpos, cigar, seq in iter:
        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        aligned = '-' * rpos
        left = 0
        for length, operation in tokens:
            length = int(length)
            if operation in 'M=X':
                aligned += seq[left:(left + length)]
                left += length
            elif operation == 'D':
                aligned += '-' * length
            elif operation in 'SI':
                left += length  # soft clip

        # pad on right
        aligned += '-'*(reflen-len(aligned))
        outfile.write('>{}\n{}\n'.format(qname, aligned))  # last entry


def encode_diffs(iter, reflen, alphabet='ACGT'):
    """
    Serialize differences of query sequences to reference
    genome, which comprise nucleotide substitutions, in-frame
    indels, and locations of missing data.
    NOTE: runs of 'N's are also represented by 'X' tokens in the CIGAR
    string.
    :param iter:  generator from minimap2()
    :param reflen:  length of reference genome
    """
    for qname, rpos, cigar, seq in iter:
        diffs = []
        missing = []
        if rpos > 0:
            # incomplete on left
            missing.append(tuple([0, rpos]))
        left = 0  # index for query

        tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
        for length, operator in tokens:
            length = int(length)
            substr = seq[left:(left + length)]
            if operator == 'X':
                # each nucleotide is a separate diff
                if 'N' in substr:
                    # for now, assume the whole substring is bs
                    missing.append(tuple([rpos, rpos+length]))
                else:
                    # assume adjacent mismatches are independent substitutions
                    for i, nt in enumerate(substr):
                        if nt in alphabet:
                            diffs.append(tuple(['~', rpos + i, nt]))
                        else:
                            # skip ambiguous base calls, like "R"
                            missing.append(tuple([rpos+i, rpos+i+1]))
                left += length
                rpos += length

            elif operator == 'S':
                # discard soft clip
                left += length

            elif operator == 'I':
                # insertion relative to reference
                diffs.append(tuple(['+', rpos, substr]))
                left += length

            elif operator == 'D':
                # deletion relative to reference
                diffs.append(tuple(['-', rpos, length]))
                rpos += length

            elif operator == '=':
                # exact match
                left += length
                rpos += length

            elif operator == 'H':
                # hard clip, do nothing
                pass
            else:
                print("ERROR: unexpected operator {}".format(operator))
                sys.exit()

        # update missing if sequence incomplete on the right
        if rpos < reflen:
            missing.append(tuple([rpos, reflen]))

        yield qname, diffs, missing


def parse_header(qname, regions, typos):
    label, accession, coldate = qname.split('|')
    if coldate.count('-') != 2:
        # collection date not full precision
        coldate = None

    try:
        _, country, _, _ = label.split('/')
    except:
        print("mangled header {}".format(label))
        return None, None, None

    if country in typos:
        country = typos[country]

    region = regions.get(country, None)
    if region is None:
        # flag devs about missing region-country entry
        print('need to update countries.json with {}'.format(country))

    # GISAID entries from China are labelled by city/province
    if region == 'China':
        country = 'China'
        region = 'Asia'

    return region, country, coldate


def parse_date(isodate):
    year, month, day = map(int, isodate.split('-'))
    return date(year, month, day)


def filter_outliers(iter, origin='2019-12-01', rate=8e-4, cutoff=0.005):
    """
    Exclude genomes that contain an excessive number of genetic differences
    from the reference, assuming that the mean number of differences increases
    linearly over time and that the variation around this mean follows a
    Poisson distribution.
    :param iter:  generator, returned by encode_diffs()
    :param origin:  str, date of root sequence in ISO format (yyyy-mm-dd)
    :param rate:  float, molecular clock rate (subs/site/yr)
    :param cutoff:  float, use 1-cutoff to compute quantile of Poisson
                    distribution
    :yield:  tuples from generator that pass filter
    """
    mu = rate*29900/365.  # convert to subs/genome/day
    t0 = parse_date(origin)
    for qname, diffs, missing in iter:
        coldate = qname.split('|')[-1]
        if coldate.count('-') != 2:
            continue
        dt = (parse_date(coldate) - t0).days
        max_diffs = poisson.ppf(q=1-cutoff, mu=mu*dt)
        ndiffs = len(diffs)
        if ndiffs > max_diffs:
            # reject genome with too many differences given date
            continue
        yield qname, diffs, missing


def collate_diffs(encoder, regions, typos, mask):
    """
    Stream output from encode_diffs to collate the incidence of each
    genetic difference by location and date, and return as a tabular
    data set.
    :param encoder:  generator, returned by encode_diffs()
    :param regions:  dict, counts keyed by region, country and collection date
    """
    res = {}
    for qname, diffs, missing in filter_outliers(encoder):
        region, country, coldate = parse_header(qname, regions, typos)
        if coldate is None:
            continue

        coldate = parse_date(coldate)
        epiweek = Week.fromdate(coldate).week
        year = coldate.year
        yeek = '{}|{:02d}'.format(year, epiweek)

        # update nested dictionaries
        if region not in res:
            res.update({region: {}})
        if country not in res[region]:
            res[region].update({country: {}})
        if yeek not in res[region][country]:
            res[region][country].update({yeek: {}})

        # iterate through genetic differences in this genome
        branch = res[region][country][yeek]  # shorthand
        for diff in diffs:
            typ, pos, alt = diff
            if typ == '~' and int(pos) in mask and alt in mask[pos]['alt']:
                # masked substitution
                continue
            if typ != '-' and 'N' in alt:
                # ignore substitutions and insertions with uncalled bases
                continue
            key = '|'.join(map(str, diff))
            if key not in branch:
                branch.update({key: 0})
            branch[key] += 1

    return res


def parse_args():
    parser = argparse.ArgumentParser("Wrapper script for minimap2")
    parser.add_argument('fasta', type=argparse.FileType('r'),
                        help="<input> path to query FASTA file")
    parser.add_argument('-o', '--outfile',
                        type=argparse.FileType('w'),
                        required=False,
                        help="<output, optional> path to write JSON output, "
                             "defaults to stdout")
    parser.add_argument('-t', '--thread', type=int, default=3, 
                        help="<option> number of threads")
    parser.add_argument('-f', '--force-headers', action='store_true',
                        help="<option> use -f to force this script to accept "
                             "headers with spaces, which will be truncated "
                             "by minimap2")
    parser.add_argument('--ref', help="<input> path to target FASTA (reference)",
                        default='data/MT291829.fa')
    parser.add_argument('--regions', default='data/countries.json',
                        help="<input> path to JSON with region to country map.")
    parser.add_argument('--typos', default='data/typos.json',
                         help="<input> path to JSON with corrections for misspelt "
                              "countries.")
    parser.add_argument('--vcf', type=str, default='data/problematic_sites_sarsCov2.vcf',
                        help="<input> path to VCF file describing problematic sites.")
    parser.add_argument('--minlen', help="<option> minimum sequence length, "
                                         "defaults to 29000nt.",
                        type=int, default=29000)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.outfile is None:
        args.outfile = sys.stdout

    # check input headers for spaces
    if not args.force_headers:
        for line in args.fasta:
            if line.startswith('>') and ' ' in line:
                print("WARNING: at least one FASTA header contains a space")
                print(line)
                print("Use '-f' to force this script to process the file")
                print("Otherwise use `sed -i 's/ /_/g' <file>` to replace all spaces in place.")
                sys.exit()

    # get length of reference
    reflen = len(convert_fasta(open(args.ref))[0][1])
    mm2 = minimap2(args.fasta.name, ref=args.ref, nthread=args.thread,
                   minlen=args.minlen)

    # load JSON files
    with open(args.regions) as fp:
        regions = json.load(fp=fp)
    with open(args.typos) as fp:
        typos = json.load(fp=fp)

    # load problematic sites from VCF
    mask = load_filter(args.vcf)

    encoder = encode_diffs(mm2, reflen=reflen)
    res = collate_diffs(encoder, regions, typos, mask)
    json.dump(res, args.outfile, sort_keys=True, indent=2, separators=(',', ': '))
