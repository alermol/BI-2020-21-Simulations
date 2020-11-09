#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import argparse
# import random
# import sys

import numpy as np
import tqdm
from Bio.Seq import Seq

import Probe
import Target


def calculate_proportion(fragment_info):
    total = fragment_info.end - fragment_info.start + 1
    try:
        vector_len = (fragment_info.adj_str["vector"][1] -
                      fragment_info.adj_str["vector"][0] + 1)
    except KeyError:
        return 1.0, 0.0
    else:
        insert_len = total - vector_len
        return (insert_len / total), (vector_len / total)


def do_simulation(simulation_number: int,
                  target_introns_number: int,
                  target_length: int = 2000,
                  targer_gc_content: float = 0.32,
                  number_of_mutated_targets: int = 1,
                  exon_hdist: int = 5,
                  inton_hdist: int = 20,
                  probe_length: int = 1000,
                  hybridization_threshold: float = 0.8,
                  probe_fragment_min_len: int = 200,
                  probe_fragment_max_len: int = 600):
    target = Target.Target(length=target_length,
                           gc_content=targer_gc_content)
    real_target = target.generate_exon_intron_target(target_introns_number)
    mt_keys = [f"mt{i}" for i in range(1, number_of_mutated_targets + 1)]
    real_target_dict = {"rt": [real_target.seq,
                               str(Seq(real_target.seq).reverse_complement())]}
    mutated_targets = {}
    for key in mt_keys:
        forward = target.mutate_target(target=real_target,
                                       exon_hdistance=exon_hdist,
                                       intron_hdistance=inton_hdist).seq
        reverse_complement = str(Seq(forward).reverse_complement())
        mutated_targets[key] = [forward, reverse_complement]
    targets_pool = {**real_target_dict, **mutated_targets}
    # generate probe from real target
    probe = Probe.Probe(target=real_target,
                        lenght=probe_length)
    probe_seq = probe.generate_probe()
    while not any(["intron" in i for i in list(probe_seq.str.keys())]):
        probe_seq = probe.generate_probe()
    # do hybridization probabilities for all fragments in pool
    specificity = {key: None for key in targets_pool.keys()}
    # for tp, target in targets_pool.items():
    true_negative = false_positive = 0
    for i in tqdm.tqdm(range(simulation_number)):
        # generate random fragment from probe from real target
        fragment_size = np.random.randint(probe_fragment_min_len,
                                          high=probe_fragment_max_len + 1)
        fragment_start = np.random.randint(0, (len(probe_seq.seq) -
                                               fragment_size))
        fragment_end = fragment_start + fragment_size - 1
        fragment_info = target.get_fragment_info(probe_seq,
                                                 fragment_start,
                                                 fragment_end)
        fragment = probe_seq.seq[fragment_start:fragment_end]
        fragment = fragment[::-1] if np.random.randint(0, 2) == 1 else fragment
        # pick target
        target_name = np.random.choice(np.asarray(list(targets_pool.keys())))
        target_seq = targets_pool[target_name]
        strand = np.random.randint(0, 2)  # 0 - forward, 1 - reverse complement
        upped_start = len(target_seq[strand]) - len(fragment)
        start_position = np.random.randint(0, upped_start)
        end_position = start_position + len(fragment)
        sequence_in_target = target_seq[strand][start_position:end_position]
        sim = sum([1 for i in zip(fragment, sequence_in_target)
                   if i[0] == i[1]]) / len(fragment)
        # true negative
        if (("rt" in target_name and
             sim < hybridization_threshold and
             (calculate_proportion(fragment_info)[1] >=
                hybridization_threshold)) or
            ("mt" in target_name and
             sim < hybridization_threshold and
             (calculate_proportion(fragment_info)[0] >=
                 hybridization_threshold))):
            true_negative += 1
        elif ("rt" in target_name and
              sim < hybridization_threshold and
              (calculate_proportion(fragment_info)[0] >=
                  hybridization_threshold)):
            false_positive += 1
    # collect data
    specificity = true_negative / (true_negative + false_positive)
    introns_in_probe = len([1 for i in probe_seq.str.keys()
                            if "intron" in i])
    introns_lengths = [v[0] - v[1] + 1 for k, v in probe_seq.str.items()
                       if "intron" in k]
    mean_inton_length_in_probe = abs(sum(introns_lengths) // introns_in_probe)
    with open("Simulation.log", "a") as log:
        log.writelines(f"{target_introns_number}\t{target_length}\t"
                       f"{targer_gc_content}\t{number_of_mutated_targets}\t"
                       f"{exon_hdist}\t{inton_hdist}\t{probe_length}\t"
                       f"{hybridization_threshold}\t{probe_fragment_min_len}\t"
                       f"{probe_fragment_max_len}\t{introns_in_probe}\t"
                       f"{mean_inton_length_in_probe}\t{specificity}\n")


if __name__ == "__main__":
    with open("Simulation.log", "w") as log:
        log.writelines("target_introns_number\ttarget_length\t"
                       "targer_gc_content\tnumber_of_mutated_targets\t"
                       "exon_hdist\tinton_hdist\tprobe_length\t"
                       "hybridization_threshold\tprobe_fragment_min_len\t"
                       "probe_fragment_max_len\tintrons_in_probe\t"
                       "mean_inton_length_in_probe\tspecificity\n")
    for i in range(500):
        try:
            do_simulation(simulation_number=100000,
                          target_introns_number=np.random.randint(1, 11),
                          target_length=np.random.randint(1000, 10000),
                          targer_gc_content=np.random.uniform(
                              low=0.2, high=0.8),
                          number_of_mutated_targets=np.random.randint(1, 11),
                          exon_hdist=np.random.randint(1, 10),
                          inton_hdist=np.random.randint(10, 100),
                          probe_length=np.random.randint(1000, 10000),
                          hybridization_threshold=0.8,
                          probe_fragment_min_len=200,
                          probe_fragment_max_len=600)
        except Exception:
            pass
