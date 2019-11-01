from pomegranate import HiddenMarkovModel


def intron_counter(seq):
    icounts = 0
    for el in seq:
        if el == 'in0 spacer0' or el == 'in1 spacer0' or el == 'in2 spacer0':
            icounts += 1
    return icounts


def find_gene_cut_index(path_names, last_elements):
    return [i for i, e in enumerate(path_names) if e in last_elements]


def find_intercoding_region(begins, ends, seq):
    cuts = []
    subseqs = []
    print(begins)
    print(ends)
    for i, begin_point in enumerate(begins):
        if i < len(ends) and i < len(begins):
            if i == 0:
                cuts.append((0, begin_point))
            else:
                cuts.append((ends[i - 1], begin_point))
    print(cuts)
    for cut in cuts:
        subseqs.append(seq[cut[0]: cut[1]])
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
    ends = find_gene_cut_index(path_names, ['stop zone taa9', 'stop zone tag9', 'stop zone tga9'])

    ext_subseq = find_intercoding_region(starts, ends, seq)

    for subs in ext_subseq:
        print('pass')
        path = predict_path(promoter_utr_model, subs)
        print([(p, subs[i - 1]) for i, p in enumerate(path) if i - 1 < len(subs)], 'h')

    #print(intron_counter(path_names))


def predict_all(seq, string):
    path_names = predict_path(coding_model, seq)
    print(path_names)
    # print([(string[i + 1], name, i - len(path_names) + 1) for i, name in enumerate(path_names) if i + 1 < len(string)])
    return path_names