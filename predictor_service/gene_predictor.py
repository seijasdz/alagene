from pomegranate import HiddenMarkovModel
import numpy
from converter_to import converter_to



def intron_counter(seq):
    icounts = 0
    for el in seq:
        if el == 'in0 spacer0' or el == 'in1 spacer0' or el == 'in2 spacer0':
            icounts += 1
    return icounts


def find_gene_cut_index(path_names, valid_elements):
    return [i for i, e in enumerate(path_names) if e in valid_elements]


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

    cds_starts = find_gene_cut_index(path_names, ['start zone7'])
    gene_end = find_gene_cut_index(path_names, ['post_poly_spacer14'])


    ext_subseq = find_intercoding_region(cds_starts, gene_end, seq)

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
    return cuts

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
    return cuts


def get_bindings(seq):
    valid_bindings = ['inr2', 'no inr2', 'tata3', '']
    bindings = []
    for i, part in enumerate(seq):
        if part in valid_bindings:
            bindings.append(i)
    return bindings


def get_exons(seq):
    in_cds = False
    cuts = []
    start = -1
    start_s = ''
    start_zones = ['utr exon', 'start zone7', 'acceptor014', 'acceptor114', 'acceptor214']
    stop_zones = ['donorx00', 'donor03', 'donor14', 'donor25', 'poly a zone 0']
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
    return cuts


def get_zones(seq):
    result = {
        'genes': []
    }
    genes = get_genes(seq)
    coding_sequences = get_cds(seq)
    exons = get_exons(seq)
    bindings = get_bindings(seq)
    for gene in genes:
        new_gene = {
            'binding': 0,
            'ss': gene,
            'exon': [],
            'cds': [],
        }

        for binding in bindings:
            if gene[0] > binding > new_gene['binding']:
                new_gene['binding'] = binding
        for exon in exons:
            if exon[0] >= gene[0] and exon[1] <= gene[1]:
                new_gene['exon'].append(exon)
        for cds in coding_sequences:
            if cds[0] >= gene[0] and cds[1] <= gene[1]:
                new_gene['cds'].append(cds)
        result['genes'].append(new_gene)
    print(result)
    return result


def predict_all(string):
    seq = numpy.array(converter_to(list(string), 2), numpy.unicode_)

    path_names = predict_path(coding_model, seq)

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

    zones = get_zones(complete_seq)

    return zones

