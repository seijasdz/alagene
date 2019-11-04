
def converter_to(sequence, order=1):
    new_list = []
    for index, element in enumerate(sequence):
        if index > order - 1:
            emission_name = element + '|'
            for x in range(order, 0, -1):
                emission_name += sequence[index - x]
            new_list.append(emission_name)
    return new_list

