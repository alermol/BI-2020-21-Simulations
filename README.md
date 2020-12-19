# Tyramide-FISH hybridization simulation

## Aim of the project

To simulate relation of hybridization events number from various parameters
of Tyramide-FISH exepriment and to find intriguing patterns.

## Methods

Tyramide-FISH is mostly used for visualization of small targets – genes or
markers. In most cases genes are included in gene families and it is tricky
to design a probe to visualize only one member of a multigenic family. The
Tyramide-FISH method has shown that including intron in a probe allows us to
visualize a specific gene from a multigenic family
([Romanov *et al*., 2015](https://pubmed.ncbi.nlm.nih.gov/26158384/)).
In our work we tried to figure out the influence of different introns on
hybridization specificity.

Tyramide-FISH method specificity of hybridization is defined by a stringency
that limits the percentage of matches and mismatches between probe and target
nucleic acid that are allowed to occur without the double helix hybrid falling
apart. The most commonly using stringency is 80%. At this stringency a hybrid
with 80% bases or more along the probe-target hybrid being complementary and
20% or less mismatched, will remain stable. Hybrid molecules with 80% or less
homology do not form or dissociate immediately.

For hybridization outcomes simulation Monte-Carlo approach is used.

Before generation script generate random target and random probe from
this target. On each iteration set of probe fragments of different length
is picked and a fraction of complementary bases between fragment and target 
is caclulated. If at least one fragment from set have equal or greater value
(fraction) threshold than iteration is successful.

## Requirements

Script writed on Python3, list of all required Python packages in
`requirements.txt` file and can be installed in virtual environment using pip

## Simulation output

As a main output of all iterations the function in script returns
a dictionary with number of hybridization events of four types: real target
hybridization only, mutate target hybridization only, both targets
hybridization and no hybridization. You can use this output as you want:
save in table (install `pandas` first), draw plots using `pyplot` and
`matplotlib` etc.

### Simulation results representation examples

#### Table form

| IN  | R   | M   | N   | B   |
| --- | --- | --- | --- | --- |
| 0   | 146 | 138 | 499 | 0   |
| 1   | 146 | 85  | 499 | 0   |
| 2   | 146 | 141 | 499 | 0   |
| 3   | 168 | 137 | 499 | 0   |
| 4   | 146 | 138 | 499 | 0   |
| 5   | 146 | 141 | 499 | 0   |
| 0   | 146 | 85  | 499 | 0   |
| 1   | 146 | 103 | 499 | 0   |
| ... | ... | ... | ... | ... |

>IN - introns number
>
>R - real target only hybridization events number
>
>M - mutated target only hybridization events number
>
>N - no hybridization events number
>
>B - both target hybridization events number

Each line is one simulation with certain number of iteration

#### Plot form

As a result of a simulation a plot like this can be created.

![example_plot](example_result/example_plot_thumb.jpeg)

## Results

In our simulation it is seen that hybridization events number with mutated
target is less than hybridization events number with real target even after
including in probe a singe intron fragment. Based on standart deviation
(whiskers on plot) we can not to confirm presence of statistically significant
difference in successful hybridization events number with real and
mutated target for intron-contained probe. This simulation will be improved
to achieve result that better match to reality.

![result_plot](results/result_plot_thumb.jpeg)

[Plot table](results/result_table.tsv)

## References

Romanov, D., Divashuk, M., Havey, M. J., & Khrustaleva, L. (2015).
Tyramide-FISH mapping of single genes for development of an integrated
recombination and cytogenetic map of chromosome 5 of Allium cepa.
Genome, 58(3), 111-119.
[[PubMed](https://pubmed.ncbi.nlm.nih.gov/26158384/)]
[[DOI](https://www.doi.org/10.1139/gen-2015-0019)]
