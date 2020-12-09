from Bio.Seq import Seq

import numpy as np
import matplotlib.pyplot as plt
import itertools


def mutate_target(sample, intron_divergency, introns_number):
    mutate_target = list(sample)
    split_positions = np.random.randint(low=0, high=len(sample),
                                        size=introns_number * 2)
    split_positions.sort()
    split_positions = np.reshape(split_positions, (2, 2))

    for i in split_positions:
        changing_letters = int((i[1] - i[0] + 1) * intron_divergency)
        pos_to_replace = np.random.randint(i[0], i[1] + 1,
                                           size=changing_letters)

        for letter in pos_to_replace:
            possible_replacement = [l for l in ('A', 'T', 'G', 'C')
                                    if l != sample[letter]]
            mutate_target[letter] = np.random.choice(possible_replacement,
                                                     size=1)[0]
    return ''.join(mutate_target)


def calculate_complementarity(seq1: str, seq2: str):
    complementary_pairs = list(zip("ATGC", "TACG"))
    total_len = len(seq1)
    complementary = 0

    for l in zip(seq1, seq2):
        if l in complementary_pairs:
            complementary += 1
    return complementary / total_len


def main(generation_number,
         target_length,
         target_gc,
         target_int_num,
         intron_divergency,
         vector_length,
         insert_length,
         min_fragment_length,
         max_fragment_length,
         complemenarity_thrs):
    target = np.random.choice(('A', 'T', 'G', 'C'),
                              size=target_length,
                              p=[(1 - target_gc) / 2,
                                 (1 - target_gc) / 2,
                                 target_gc / 2,
                                 target_gc / 2])
    target = ''.join(target)

    mutated_target = mutate_target(target, intron_divergency, target_int_num)
    target_pool = [target, mutated_target]

    vector = np.random.choice(('A', 'T', 'G', 'C'),
                              size=vector_length)
    vector = ''.join(vector)
    insert_start = int(np.random.randint(0, target_length - insert_length + 1))
    insert_end = int(insert_start + insert_length)
    probe = target[insert_start:insert_end]
    probe = vector + probe + vector[:max_fragment_length]

    result = np.zeros(generation_number, dtype=np.float32)
    for n in range(generation_number):
        successful_hybrisization = [None, None]
        for idx, t in enumerate(target_pool):
            fragments = []
            while sum(list(map(len, fragments))) <= target_length:
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
                'N' for _ in range(target_length - sum(list(map(len,
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
            result[n] = 1.0  # real target only
        elif (successful_hybrisization[0] == False and
              successful_hybrisization[1] == True):
            result[n] = -1.0  # mutated target only
        elif (successful_hybrisization[0] == False and
              successful_hybrisization[1] == False):
            result[n] = 0.0  # no hybridization
        else:
            result[n] = -0.5 # both targets
    return (np.size(result[result == 1.0]),
            np.size(result[result == -1.0]),
            np.size(result[result == 0.0]),
            np.size(result[result == -0.5]))


if __name__ == "__main__":
    np.random.seed(5671)
    data = main(generation_number=50000,
                target_length=3000,
                target_gc=0.5,
                target_int_num=2,
                intron_divergency=0.5,  # fraction of divergent bases in introns
                vector_length=3000,
                insert_length=1500,
                min_fragment_length=200,
                max_fragment_length=600,
                complemenarity_thrs=0.8)
    print(
f'''
Experiment result:

Hybridization with real target only: {data[0]}
Hybridization with mutated target only: {data[1]}
Hybridization with both targets: {data[2]}
No hybridization with both targets: {data[3]}
'''
    )
