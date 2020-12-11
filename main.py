'''
Author: Aleksey Ermolaev

Script for simulation of hybridization of DNA-probe with DNA-target
in Tyramide-FISH method
'''

from Bio.Seq import Seq

import numpy as np
import matplotlib.pyplot as plt
import itertools


def mutate_hyb_site(site_seq, introns_divergency, introns_number):
    '''Return sequence with changed fraction of letters in introns (mutated)

    Keyword arguments:
    site_seq -- sequence to mutate (str)
    introns_divergency -- fraction of letters to change in each intron (float)
    introns_number -- number of areas to intron location (int)
    '''
    mutated_site_seq = list(site_seq)
    split_positions = sorted(np.random.randint(low=0, high=len(site_seq),
                                               size=introns_number * 2))
    # pairs of values - introns boundaries
    split_positions = np.reshape(split_positions, (introns_number, 2))

    for i in split_positions:
        changing_letters = int((i[1] - i[0] + 1) * introns_divergency)
        pos_to_replace = np.random.randint(i[0], i[1] + 1,
                                           size=changing_letters)

        for letter in pos_to_replace:
            possible_replacement = {'A', 'T', 'G', 'C'}
            possible_replacement.remove(site_seq[letter])
            possible_replacement = list(possible_replacement)
            mutated_site_seq[letter] = np.random.choice(possible_replacement,
                                                        size=1)[0]
    return ''.join(mutated_site_seq)


def calculate_complementarity(seq1, seq2):
    '''Return percent of complementary pairs in two sequences

    Keyword arguments:
    seq1 -- first sequence (str)
    seq2 -- second sequence (str)

    Both first and second sequences must be an equal length
    '''
    complementary_pairs = list(zip("ATGC", "TACG"))
    total_len = len(seq1)
    complementary = 0

    for l in zip(seq1, seq2):
        if l in complementary_pairs:
            complementary += 1
    return complementary / total_len


def main(iterations_number,
         site_length,
         gc_content,
         site_intron_number,
         intron_divergency,
         vector_length,
         min_fragment_length,
         max_fragment_length,
         complemenarity_thrs):
    '''Does simualtion iterations and return the result

    Keyword arguments:
    iteration_number -- number of iterations to select fragments (int)
    site_length -- length of hybridization site (int)
    gc_content -- GC-content of hybridization site (float)
    site_intron_number -- number of introns in hybridization site (int)
    intron_divergency -- fraction of differ letters in introns (float)
    vector_length -- length of supplementary sequence for probe (int)
    min_fragment_length -- minimal length of probe fragment (int)
    max_fragment_length -- maximal length of probe fragment (int)
    complemenarity_thrs -- threshold of complementarity for stable hybrid
                           formation (float)
    '''
    prob = [(1 - gc_content) / 2] * 2 + [gc_content / 2] * 2
    hybridization_site = np.random.choice(('A', 'T', 'G', 'C'),
                                          size=site_length,
                                          p=prob)
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
    for n in range(iterations_number):
        successful_hybridization = [None, None]
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
                successful_hybridization[idx] = True
            else:
                successful_hybridization[idx] = False

        if (successful_hybridization[0] and
                not successful_hybridization[1]):
            result['R'] += 1  # real target only
        elif (not successful_hybridization[0] and
              successful_hybridization[1]):
            result['M'] += 1  # mutated target only
        elif (not successful_hybridization[0] and
              not successful_hybridization[1]):
            result['N'] += 1  # no hybridization
        else:
            result['B'] += 1  # both targets
    return result


if __name__ == "__main__":
    np.random.seed(5671)
    data = main(iterations_number=10,
                site_length=1100,
                gc_content=0.32,
                site_intron_number=2,
                intron_divergency=0.5,
                vector_length=3000,
                min_fragment_length=100,
                max_fragment_length=1100,
                complemenarity_thrs=0.8)
