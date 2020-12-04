from Bio import SeqIO
import numpy as np
import re


def read_fasta_alignment(filename, max_gap_fraction=1.0, header_regex=None):
    inds = []
    seqs = []
    fseqlen = 0
    for i, record in enumerate(SeqIO.parse(filename, "fasta")):
        i += 1
        ngaps = 0
        seq = record.seq

        if i == 1:
            fseqlen = len(seq)
            for j in range(len(seq)):
                c = seq[j]
                if c != '.' and c == c.upper():
                    inds.append(j)
                    ngaps = ngaps + 1 if c == '-' else ngaps
        else:
            if len(seq) != len(inds):
                raise Exception("Inputs are not aligned")
            for j in range(len(seq)):
                c = seq[j]
                if c != '.' and c == c.upper():
                    if inds[j] != j:
                        raise Exception("Inconsistent inputs")
                    ngaps = ngaps + 1 if c == '-' else ngaps
            if fseqlen != len(seq):
                raise Exception("Inconsistent inputs")

        if ngaps / fseqlen <= max_gap_fraction:
            seqs.append(i)

    if len(seqs) is None:
        raise Exception(f"out of {i} sequences, none passed the filter (max_gap_fraction={max_gap_fraction})")

    # Z = [[None for i in range(len(seqs))] for j in range(fseqlen)]
    Z = np.empty((fseqlen, len(seqs)))
    header = []
    sequence = []
    alphabet = [1, 21, 2, 3, 4, 5, 6, 7, 8, 21, 9, 10, 11, 12, 21, 13, 14, 15, 16, 17, 21, 18, 19, 21, 20]

    for seqid, record in enumerate(SeqIO.parse(filename, "fasta")):
        header.append(record.name)
        sequence.append(record.seq)
        if seqs[-1] < seqid + 1:
            break
        if seqs[seqid] == seqid:
            continue
        for i in range(fseqlen):
            c = record.seq[inds[i]]
            Z[i][seqid] = letter2num(alphabet, c)
    assert seqid == len(seqs) - 1

    return (Z.shape[0], Z.shape[1], np.max(Z), Z.T, sequence, header, _, _, _)


def letter2num(alphabet, c):
    i = ord(c) - ord('A')
    if 0 <= i <= 24:
        return alphabet[i]
    return 21


def specname(s, header_regex, captureinds):
    if header_regex is not None:
        res = re.search(header_regex, s)
        if res:
            captures = res.groups()
            if len(captures) < 2:
                raise Exception(
                    f"invalid header regex: should always return at least 2 captured groups if it matches; has returned: {len(captures)}")
            uniprot_id, spec_name = captures[captureinds[0]], captures[captureinds[1]]

            if not isinstance(uniprot_id, str):
                raise Exception(f"the capture group for `id` did not match the spec string: {s}")
            if not isinstance(spec_name, str):
                raise Exception(f"the capture group for `species` did not match the spec string: {s}")
        else:
            raise Exception(f"unrecognized spec string: {s}")
    else:
        regex_uniprot = r"^(?:(?:(?:[^|/_]+?\|){2})|(?:[^|/_]+?/))([^|/_]+?)_([^|/_\s]+)"
        regex_oldskrr = r"\[(.*?)\]"
        regex_joined = r"^(.*?)::(.*?)/(.*)$"

        # standart format
        res_uniprot = re.search(regex_uniprot, s)
        res_oldskrr = re.search(regex_oldskrr, s)
        res_joined = re.search(regex_joined, s)
        if res_uniprot:
            uniprot_id, spec_name = res_uniprot.groups()
        elif res_oldskrr:
            spec_name = res_oldskrr.groups()[0]
            uniprot_id = "000000"
        elif res_joined:
            spec_name = res_joined.groups()[2]
            uniprot_id = "000000"
        else:
            raise Exception(f"unrecognized spec string: {s}")
    return uniprot_id, spec_name


def compute_spec(header, header_regex=None):
    spec_name = []
    uniprot_id = []

    captureinds = (0, 1)

    if header_regex is not None:
        cnames = capture_names(header_regex)
        if "id" in cnames and "species" in cnames:
            captureinds = (cnames.index("id"), cnames.index("species"))
        elif "id" not in cnames and "species" not in cnames:
            pass
        else:
            raise Exception(
                "only one of `id` and `species` capture names found in `header_regex`; either none or both should be used")

    for i in range(len(header)):
        uid, sname = specname(header[i], header_regex, captureinds)
        uniprot_id.append(uid)
        spec_name.append(sname)

    specunique = set(spec_name)
    sdict = {s: i for i, s in enumerate(specunique)}
    specid = [sdict[sn] for sn in spec_name]
    return specid, spec_name, uniprot_id


def capture_names(header_regex):
    # reg = r"\<.*?\>"
    reg = r"\<(.*?)\>|(\/\(\[\^\/\])"
    res = re.findall(reg, header_regex)
    names = [raw_name[0] for raw_name in res if "/" not in raw_name[0] ]
    return names


if __name__ == '__main__':
    # alphabet = [1, 21, 2, 3, 4, 5, 6, 7, 8, 21, 9, 10, 11, 12, 21, 13, 14, 15, 16, 17, 21, 18, 19, 21, 20]
    # print(letter2num(alphabet, 'C'))
    # read_fasta_alignment("../data/X1.fasta")
    header = ["SOMELOCATION/SOMESPECIES/OTHERINFO"]
    answer = ([0], ['SOMESPECIES'], ['SOMELOCATION'])

    assert compute_spec(header, r"^([^/]+)/([^/]+)/") == answer
    assert compute_spec(header, r"^([^/]+)/([^/]+)/(.*)") == answer
    assert compute_spec(header, r"^(?P<id>[^/]+)/(?P<species>[^/]+)/") == answer
    assert compute_spec(header, r"^(?P<id>[^/]+)/(?P<species>[^/]+)/") == answer
    assert compute_spec(header, r"^(?P<id>[^/]+)/(?P<species>[^/]+)/(?P<rest>.*)") == answer

    header = ["SOMESPECIES/OTHERINFO/SOMELOCATION"]
    assert compute_spec(header, r"^(?P<species>[^/]+)/([^/]+)/(?P<id>.+)") == answer
    # header_regex = r"^(?P<id>[^\/]+)\/(?P<species>[^\/]+)\/(?P<rest>.*)"
    # print(compute_spec(header, header_regex))
