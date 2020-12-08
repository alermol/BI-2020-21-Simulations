from Bio.Seq import Seq

import numpy as np
import matplotlib.pyplot as plt
import itertools


def mutate_target(sample, hdistance):
    mutated_sequence = list(sample)
    positions_to_mutate = np.random.randint(low=0,
                                            high=len(mutated_sequence) - 1,
                                            size=hdistance)
    for position in positions_to_mutate:
        original_letter = mutated_sequence[position]
        possible_replacement = [l for l in ('A', 'T', 'G', 'C')
                                if l != original_letter]
        mutated_sequence[position] = np.random.choice(possible_replacement,
                                                      size=1)[0]
    return ''.join(mutated_sequence)


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
         mtarget_hdist,
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

    mutated_target = mutate_target(target, mtarget_hdist)
    target_pool = [target, mutated_target]

    vector = np.random.choice(('A', 'T', 'G', 'C'),
                              size=vector_length)
    vector = ''.join(vector)
    insert_start = int(np.random.randint(0, target_length - insert_length + 1))
    insert_end = int(insert_start + insert_length)
    probe = target[insert_start:insert_end]
    probe = vector + probe + vector[:max_fragment_length]

    result = []  # array to return any data from simulation function
    for n in range(generation_number):
        for t in target_pool:
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
    return result


if __name__ == "__main__":
    np.random.seed(5671)
    data = main(generation_number=100,
                target_length=3000,
                target_gc=0.5,
                mtarget_hdist=200,
                vector_length=3000,
                insert_length=1500,
                min_fragment_length=200,
                max_fragment_length=600,
                complemenarity_thrs=0.8)
