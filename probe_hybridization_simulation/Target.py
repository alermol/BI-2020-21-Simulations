import itertools
from collections import namedtuple

import numpy as np


class Target:
    def __init__(self,
                 length: int,
                 gc_content: float):
        self.length = length
        self.gc_content = gc_content
        self.ALPHABET = ["A", "T", "G", "C"]
        assert self.length >= 700, "Target size too small"
        assert 0 <= self.gc_content <= 1, "Wrong GC-content value"


    def __generate_output(self,
                          sequence: str,
                          sequence_type: str,
                          structure: dict,
                          mutated: bool):
        output = namedtuple("Target", "seq type str mutated")
        return output(sequence, sequence_type, structure, mutated)

    def generate_exon_target(self):
        structure = {}
        gc_probability = self.gc_content / 2
        at_probability = (1 - self.gc_content) / 2
        nucl_probability = ([at_probability] * 2) + ([gc_probability] * 2)
        sequence = np.random.choice(self.ALPHABET,
                                    size=self.length,
                                    p=nucl_probability)
        structure["exon1"] = [1, len(sequence)]
        output = self.__generate_output(sequence="".join(sequence),
                                        sequence_type="intonless",
                                        structure=structure,
                                        mutated=False)
        return output


    def generate_exon_intron_target(self, intron_number: int = 1):
        assert intron_number > 0, ("Introns number can not be 0\n"
                                   "For generate intronless probe use "
                                   "generate_exon_target() function")
        structure = {}
        gc_probability = self.gc_content / 2
        at_probability = (1 - self.gc_content) / 2
        nucl_probability = ([at_probability] * 2) + ([gc_probability] * 2)
        sequence = np.random.choice(self.ALPHABET,
                                    size=self.length,
                                    p=nucl_probability)
        exon_boundary = np.random.random_integers(3, high=(self.length - 2),
                                                  size=intron_number * 2)
        boundaries = sorted(exon_boundary)
        # add first and last positions number in boundaries map
        boundaries.insert(0, 1)
        boundaries.append(self.length)
        # generate numbers for intons and exons for structure
        numbers_structure = [[i] * 2 for i in range(1, len(boundaries) + 1)]
        del numbers_structure[-1]
        del numbers_structure[-1][-1]
        numbers_structure = list(itertools.chain(*numbers_structure))
        # generate exon-inron structure names for structures
        exon_intron_structure = [("exon", "intron") for _ in
                                 range(len(boundaries) // 2)]
        exon_intron_structure = list(itertools.chain(*exon_intron_structure))
        exon_intron_structure += ["exon"]
        exon_intron_structure = [f"{seq_type}{number}" for seq_type, number in
                                 zip(exon_intron_structure, numbers_structure)]
        # generate bounaries pairs for each exon and intron
        boundaries_coord = [(boundaries[i], boundaries[i + 1]) for i in
                            range(0, len(boundaries) - 1)]
        boundaries_coord = [(value[0] + 1, value[1] - 1) if index % 2 == 1
                             else (value[0], value[1]) for index, value in
                             enumerate(boundaries_coord)]
        # generate final structure
        structure = {name: coord for name, coord in zip(exon_intron_structure,
                                                        boundaries_coord)}
        output = self.__generate_output(sequence="".join(sequence),
                                        sequence_type="contain_introns",
                                        structure=structure,
                                        mutated=False)
        return output

    @staticmethod
    def get_position_info(self, target, position: int):
        assert 1 <= position <= len(target.seq), "Position out of sequence"
        nucleotide = target.seq[position - 1]
        area_annotation = {}
        for annotation, boundary in target.str.items():
            if boundary[1] >= position >= boundary[0]:
                area_annotation[annotation] = boundary
                break
            else:
                continue
        output = namedtuple("Position_info", "pos nucl annot")
        return output(pos=position, nucl=nucleotide, annot=area_annotation)


    @staticmethod
    def get_fragment_info(target, start: int, end: int):
        assert 1 <= start <= len(target.seq), "Position out of sequence"
        assert 1 <= end <= len(target.seq), "Position out of sequence"
        assert start <= end, "Start position greater or equal than end position"
        # find annotation of area for start position
        start_area = ""
        for annotation, boundary in target.str.items():
            if boundary[1] >= start >= boundary[0]:
                start_area += annotation
                break
            else:
                continue
        # find annotation of area for end position
        end_area = ""
        for annotation, boundary in target.str.items():
            if boundary[1] >= end >= boundary[0]:
                end_area += annotation
                break
            else:
                continue
        # create area annotation for fragment
        start_area_index = list(target.str.keys()).index(start_area)
        end_area_index = list(target.str.keys()).index(end_area) + 1
        selected_areas = list(target.str.keys())[start_area_index:end_area_index]
        area_annotation = {area:target.str[area] for area in selected_areas}
        # generate output
        output = namedtuple("Fragment_info",
                            "seq start end length inum enum str")
        return output(seq=target.seq[start:end],
                      start=start,
                      end=end,
                      length=len(target.seq[start:end]),
                      inum=sum(["intron" in i for i in selected_areas]),
                      enum=sum(["exon" in i for i in selected_areas]),
                      str=area_annotation)


    def mutate_target(self,
                      target,
                      exon_hdistance: int,
                      intron_hdistance: int):
        """
        Hamming distance will pe applied to each intron and exon.
        Summary hamming distance for mutated probe can be calculated as:

        (exon_hdistance * exon_number) + (intron_hdistance * intron_number)
        """
        assert (exon_hdistance > 0 or
                intron_hdistance > 0), "Hamming distance must be >0"
        exon_number = sum(["exon" in i for i in target.str.keys()])
        intron_number = sum(["intron" in i for i in target.str.keys()])
        general_hdist = (exon_hdistance * exon_number +
                         intron_hdistance * intron_number)
        assert general_hdist <= len(target.seq), "Hamming distance is too big"
        original_seq = list(target.seq)
        for area, coord in target.str.items():
            if "exon" in area:
                mutated_positions = np.random.random_integers(
                    coord[0],
                    high=coord[1],
                    size=exon_hdistance)
            else:
                mutated_positions = np.random.random_integers(
                    coord[0],
                    high=coord[1],
                    size=intron_hdistance)
            for pos in mutated_positions:
                original_letter = target.seq[pos - 1]
                possible_replace = [i for i in self.ALPHABET
                                    if i != original_letter]
                original_seq[pos - 1] = np.random.choice(possible_replace,
                                                     size=1)[0]
        return self.__generate_output(sequence="".join(original_seq),
                                      sequence_type=target.type,
                                      structure=target.str,
                                      mutated=True)
