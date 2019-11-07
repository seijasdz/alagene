from pomegranate import HiddenMarkovModel
import numpy
from converter_to import converter_to


def intron_counter(seq):
    icounts = 0
    for el in seq:
        if el == 'in0 spacer0' or el == 'in1 spacer0' or el == 'in2 spacer0':
            icounts += 1
    return icounts


def find_gene_cut_index(path_names, last_elements):
    return [i for i, e in enumerate(path_names) if e in last_elements]


def find_intercoding_region(begins, ends, seq):
    subseqs = []
    for i, b in enumerate(begins):
        if i == 0:
            before = 0
        else:
            before = ends[i - 1] + 1
        subseqs.append(seq[before: b])
    return subseqs


def predict_path(model, seq):
    logp, path = model.viterbi(seq)
    return [p[1].name for p in path]


with open('coding_model_base_poly.json') as base_model_file:
    coding_model_json = base_model_file.read()

with open('utr_model_base.json') as promoter_model_file:
    promoter_utr_model_json = promoter_model_file.read()

coding_model = HiddenMarkovModel.from_json(coding_model_json)
promoter_utr_model = HiddenMarkovModel.from_json(promoter_utr_model_json)


def predict_all_old(seq, string):
    path_names = predict_path(coding_model, seq)

    print([(string[i + 1], name, i - len(path_names) + 1) for i, name in enumerate(path_names) if i + 1 < len(string)])

    starts = find_gene_cut_index(path_names, ['start zone7'])
    ends = find_gene_cut_index(path_names, ['post_poly_spacer14'])

    ext_subseq = find_intercoding_region(starts, ends, seq)

    for subs in ext_subseq:
        print('pass')
        path = predict_path(promoter_utr_model, subs)
        print([(p, subs[i - 1]) for i, p in enumerate(path) if i - 1 < len(subs)], 'h')

    #print(intron_counter(path_names))


def divide_genes(begins, ends, seq):
    divided = []
    for i, b in enumerate(begins):
        if i < len(ends):
            end = ends[i] + 1
        else:
            end = -1
        divided.append(seq[b: end])
    return divided


def get_genes(seq):
    in_gene = False
    cuts = []
    start = -1
    for i, part in enumerate(seq):
        if part != 'back' and not in_gene:
            in_gene = True
            start = i
        elif part == 'back' and in_gene:
            in_gene = False
            stop = i
            cuts.append((start, stop))
        elif part != 'back' and in_gene and i == len(seq) -1:
            stop = i
            cuts.append((start, stop))
    print(cuts)

def get_cds(seq):
    in_cds = False
    cuts = []
    start = -1
    start_s = ''
    start_zones = ['start zone7', 'acceptor014', 'acceptor114', 'acceptor214']
    stop_zones = ['donor03', 'donor14', 'donor25', 'stop zone tga9', 'stop zone tag9', 'stop zone taa9']
    for i, part in enumerate(seq):
        if part in start_zones and not in_cds:
            in_cds = True
            start = i
            start_s = part
        elif part in stop_zones and in_cds:
            in_cds = False
            cuts.append((start, i))
        elif in_cds and i == len(seq) - 1:
            cuts.append((start, i))
    print(cuts)

def get_zones(seq):
    genes = get_genes(seq)
    get_cds(seq)


def predict_all(string):
    seq = numpy.array(converter_to(list(string), 2), numpy.unicode_)

    path_names = predict_path(coding_model, seq)
    print(path_names)
    starts = find_gene_cut_index(path_names, ['start zone7'])
    ends = find_gene_cut_index(path_names, ['post_poly_spacer14'])

    genes = divide_genes(starts, ends, path_names)
    utr5s = find_intercoding_region(starts, ends, seq)

    full_seqs = []
    complete_seq = []
    for i, utr in enumerate(utr5s):
        path = predict_path(promoter_utr_model, utr)
        complete = path[1:-1] + genes[i]
        # full_seqs.append(complete)
        complete_seq += complete
    get_zones(complete_seq)
    #print(len(complete_seq), len(seq))
    #print([(x, seq[i]) for i, x in enumerate(complete_seq)])
    # print([(string[i + 1], name, i - len(path_names) + 1) for i, name in enumerate(path_names) if i + 1 < len(string)])
    return path_names

