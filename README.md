# neuralAnalysis_mat

Repository to investigate activity pattenrs of electrophysiological neural data recorded by invasive brain implants in mice using Matlab. 

## Analytical methods used: 

- **Pearson's correlation**: used to compare activity patterns of two neurons (single cell correlations). Can also be used to compare population activity between two different positions within an environment (population vector analysis). Based on functions by [Barry lab](https://github.com/Barry-lab).
- **Bayesian decoding**: computational algorithm that infers the most likely underlying state (behaviour) from observed data (neural spikes) by incorporating prior knowledge and updating beliefs using **Bayesian inference** principles. Based on a tutorial by the [Van der Meer lab](https://github.com/vandermeerlab).
- **Non-parametric statistical tests**: Wilcoxon ranksum test and Kruskal-Wallis, combined with a Bonferroni-corrected multiple comparison posthoc test, were used for **hypothesis testing**. 

## Description: 

- **doPositionalShift.m**: function to calculate any potential change in position of a place field across trials and/or across rooms for each place cell. 
- **doRateRemapping.m**: function to calculate any potential change in firing rate of a place field across trials and/or across rooms for each place cell. 
- **performBayesianDecoding.m**: function to decode the position of the animal based on activity patterns of place cells within a given session. 
- **populationVectorAnalysis.m**: function to compute the correlation in population activity across spatial bins in an environment.
- **singleCellCorrelation.m**: function to compute the correlation in activity pattern across place cells in an environment. 

## Glossary:

- Place cell: neuron in the hippocampus that fires in a specific location in space.
- Place field: location in which a place cell is active.
- Ratemap: average activity (aka firing rate in spikes per second) of a neuron in each spatial bin of an environment explored by an animal. If the environment is a linear track, the ratemap will be a 1-D vector. If the environment is an open arena, the ratemap will be a 2-D matrix.
- Hippocampus: brain structure involved in memory and spatial navigation.
