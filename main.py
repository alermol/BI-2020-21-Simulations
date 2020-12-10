from Bio.Seq import Seq

import numpy as np
import matplotlib.pyplot as plt
import itertools

import sys
import tqdm


def mutate_hyb_site(site_seq, introns_divergency, introns_number):
    mutated_site_seq = list(site_seq)
    split_positions = np.random.randint(low=0, high=len(site_seq),
                                        size=introns_number * 2)
    split_positions.sort()
    split_positions = np.reshape(split_positions, (introns_number, 2))

    for i in split_positions:
        changing_letters = int((i[1] - i[0] + 1) * introns_divergency)
        pos_to_replace = np.random.randint(i[0], i[1] + 1,
                                           size=changing_letters)

        for letter in pos_to_replace:
            possible_replacement = [l for l in ('A', 'T', 'G', 'C')
                                    if l != site_seq[letter]]
            mutated_site_seq[letter] = np.random.choice(possible_replacement,
                                                        size=1)[0]
    return ''.join(mutated_site_seq)


def calculate_complementarity(seq1: str, seq2: str):
    complementary_pairs = list(zip("ATGC", "TACG"))
    total_len = len(seq1)
    complementary = 0

    for l in zip(seq1, seq2):
        if l in complementary_pairs:
            complementary += 1
    return complementary / total_len


def main(generation_number,
         site_length,
         gc_content,
         site_intron_number,
         intron_divergency,
         vector_length,
         min_fragment_length,
         max_fragment_length,
         complemenarity_thrs):
    hybridization_site = np.random.choice(('A', 'T', 'G', 'C'),
                                          size=site_length,
                                          p=[(1 - gc_content) / 2,
                                             (1 - gc_content) / 2,
                                             gc_content / 2,
                                             gc_content / 2])
    hybridization_site = ''.join(hybridization_site)

    mutated_hyb_site = mutate_hyb_site(hybridization_site,
                                       intron_divergency,
                                       site_intron_number)
    hybrid_sites = [hybridization_site, mutated_hyb_site]

    vector = np.random.choice(('A', 'T', 'G', 'C'),
                              size=vector_length)
    vector = ''.join(vector)
    probe = vector + hybridization_site + vector[:max_fragment_length]

    result = {'R': 0, 'M': 0, 'B': 0, 'N': 0}
    for n in tqdm.tqdm(range(generation_number)):
        successful_hybrisization = [None, None]
        for idx, t in enumerate(hybrid_sites):
            fragments = []
            while sum(list(map(len, fragments))) <= site_length:
                length = np.random.randint(min_fragment_length,
                                           max_fragment_length + 1)
                length = int(length)
                start = np.random.randint(low=0,
                                          high=len(probe) - length + 1,
                                          size=1)
                start = int(start)
                fragments.append(probe[start:(start + length + 1)])

            fragments = fragments[:-1]
            ns = [
                'N' for _ in range(site_length - sum(list(map(len,
                                                              fragments))))
            ]
            fragments.extend(ns)
            np.random.shuffle(fragments)
            fragments = [
                ''.join(list(i)) for _, i in itertools.groupby(fragments)
            ]

            pos_fragments = {}
            for f in range(len(fragments)):
                start = sum(map(len, fragments[:f]))
                end = start + len(fragments[f]) - 1
                key = (start, end)
                value = str(Seq(fragments[f]).complement())
                pos_fragments[key] = value

            hybridizaton_res = {}
            for pos, seq in pos_fragments.items():
                target_fragment = t[pos[0]:(pos[1] + 1)]
                compl = calculate_complementarity(seq, target_fragment)
                if compl > 0:
                    hybridizaton_res[pos] = compl

            if any([True for i in hybridizaton_res.values()
                    if i >= complemenarity_thrs]):
                successful_hybrisization[idx] = True
            else:
                successful_hybrisization[idx] = False

        if (successful_hybrisization[0] == True and
                successful_hybrisization[1] == False):
            result['R'] += 1  # real target only
        elif (successful_hybrisization[0] == False and
              successful_hybrisization[1] == True):
            result['M'] += 1  # mutated target only
        elif (successful_hybrisization[0] == False and
              successful_hybrisization[1] == False):
            result['N'] += 1  # no hybridization
        else:
            result['B'] += 1  # both targets
    return result


if __name__ == "__main__":
    # np.random.seed(5671)

    int_numbs = np.arange(0, 6, 1, dtype=int)
    ry_array = np.zeros(int_numbs.shape)
    my_array = np.zeros(int_numbs.shape)

    internal_iterations = 50

    for i in tqdm.tqdm(int_numbs):
        ry_array_tmp = np.zeros(internal_iterations, dtype='uint8')
        my_array_tmp = np.zeros(internal_iterations, dtype='uint8')
        intr_numb = np.full(internal_iterations, i)
        for ii in tqdm.tqdm(range(internal_iterations)):
            data = main(generation_number=1000000,
                        site_length=1100,
                        gc_content=0.32,
                        site_intron_number=i,
                        intron_divergency=0.5,  # fraction of divergent bases in introns
                        vector_length=3000,
                        min_fragment_length=100,
                        max_fragment_length=1100,
                        complemenarity_thrs=0.2)
            ry_array_tmp[ii] = data['R']
            my_array_tmp[ii] = data['M']
        ry_array[i] = np.mean(ry_array_tmp)
        my_array[i] = np.mean(my_array_tmp)
        plt.errorbar(i, np.mean(ry_array_tmp), ecolor='green',
                     yerr=np.std(ry_array_tmp), elinewidth=1, capsize=3)
        plt.errorbar(i, np.mean(my_array_tmp), ecolor='blue',
                     yerr=np.std(my_array_tmp), elinewidth=1, capsize=3)

    plt.plot(int_numbs, ry_array, c='green', marker='.', linestyle='--')
    plt.plot(int_numbs, my_array, c='blue', marker='.', linestyle='--')

    plt.grid()
    plt.legend(['Успешные гибридизации с реальным таргетом',
                'Успешные гибридизации с мутантным таргетом'],
               loc='upper center', fontsize='xx-small',
               bbox_to_anchor=(0.5, 1.06), shadow=True, ncol=2)
    plt.title(('Зависимость количества событий гибридизации\n'
               'от количества интронов в пробе'), pad=17.0)
    plt.xticks(int_numbs)
    plt.xlabel('Количество интронов')
    plt.ylabel('Количество событий успешной гибридизации')
    plt.savefig('plot5.png', dpi=600)
    plt.close()

    # print(
# f'''
# Experiment result:
#
# Hybridization with real target only: {data[0]}
# Hybridization with mutated target only: {data[1]}
# Hybridization with both targets: {data[3]}
# No hybridization with both targets: {data[2]}
# '''
    # )
