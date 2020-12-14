# BI-2020-21-Simulations

## Aim of the project

To simulate relation of hybridization events number from various parameters
of Tyramide-FISH exepriment and to find intriguing patterns.

## Methods

For hybridization outcomes simulation Monte-Carlo simulation is use.

## Requirements

Script writed on Python3, list of all required Python packages in
`requirements.txt` file and can be installed in virtual environment using pip

## Simulation output

As a main output of all iterations the function in script is returns
a dictionary with number of hybridization events of four types: real target
hybridization only, mutate target hybridization only, both targets
hybridization only and no hybridization. You can use this output as you want:
save in table (install `pandas` first), draw plots using `pyplot` and
`matplotlib` etc.

## References

Article about influences of introns in probe at probe specificity

Romanov, D., Divashuk, M., Havey, M. J., & Khrustaleva, L. (2015).
Tyramide-FISH mapping of single genes for development of an integrated
recombination and cytogenetic map of chromosome 5 of Allium cepa.
Genome, 58(3), 111-119.
