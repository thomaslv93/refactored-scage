# Refactored Version of [scAge](https://github.com/alex-trapp/scAge)
Refactored version of Alex Trapp's [scAge](https://github.com/alex-trapp/scAge).

## Overview of Epigenetic Clocks

### CpG Sites
CpG sites are sites along the genome in which a cytosine is immediately followed by a guanine. In this pair, the cytosine can be methylated, meaning that a methyl group can be added to it. A CpG site can be in a methylated or unmethylated state, and there are processes which can methylate an unmethylated site and processes which can demethylate a methylated site. The state of a CpG site's methylation can affect the gene expression of the surrounding DNA sequence. The CpG sites themselves--meaning their existence and location in the genome--remain consistent across an organism's lifespan, and two organisms of the same species will have essentially identical CpG sites, just as they have essentially identical genomes. The methylation state of a CpG site can be different between two organisms, within the same organism at two different times, within the same organism in two different cells from two different tissues, or even within the same organism in two different cells in the same tissue. The overall methylation state of these CpG sites in an organism is called the methylome.

### Predicting Age
The methylome exhibits predictable changes with age. These age-related changes have been used to produce "epigenetic clocks". Trained on methylation data paired with actual "chronological" ages of samples, epigenetic clocks can take in methylation data from a handful of CpG sites and predict the age of a given sample. This age is interpreted as a "biological" age. In contrast to "chronological" age which represents the actual amount of time an organism has spent on this earth, "biological" age represents an organisms point along the underlying biolical process. If someone's biological age is below their chronological age, they are considered to have a lower rate of aging. If their biological age is above their chronological age, they are considered to have a higher rate of aging.

### Multi-Cell vs Single-Cell
Individual CpG sites can be either methylated or unmethylated. However, usually methylation data is extracted from bulk samples, meaning we get reads from many cells and can then determine what fraction of these cells are methylated at different CpG sites. This yields a number between 0 and 1. This fraction varies with some predictability along the interval `[0, 1]` as an organism ages. Most epigenetic clocks use ElasticNet regression to select a subset of CpG sites and then perform a linear transformation from their fractional methylation values to a predicted age. Clocks developed in this manner are not effective for single cell samples, in which methylation values are either zero or one. Additionally, when measuring the methylation of single cells we are not guaranteed to get reads from every CpG site. Therefore the available CpG sites change from sample to sample. Alex Trapp's scAge is able to overcome these challenges by

## Overview of scAge
1. For a large number of CpG sites, use bulk-sample training data to create a linear model from age to methylation level for each CpG site. Given an age, this should return a number between zero and one. Interpret this value as the probability that the given CpG will be methylated at the input age. One minus this value will represent the probability that the given CpG will be unmethylated at the input age.
2. For simplicity, assume that CpG sites change independently. Now given methylation values of 0 or 1 for a handful of CpG sites, one can compute the probability that these methylation values occur simultaneously at the input age by multiplying the probabilities for each CpG site computed in (1). Many CpG sites have no correlation with age, so we will filter for the most age-associated CpG sites first.
3. Compute this joint probability for each age in a range of valid ages. Then select the age which maximizes the probability of a given methylome. This age is the predicted age of the sample.

## Refactoring in Response to Reviewer Comments
- Made code more compact, modular, and functional to make tweaks and modifications easier and faster.
- Increased flexibilty of the code to allow for faster experimentation with different similarity metrics.
- Improved efficiency of algorithm by converting slow for-loops to efficient matrix operations using pandas and numpy.
- Made code more efficient by filtering 1% of rows directly from file, instead of reading entire file into memory and then filtering rows.
- Wrote code to efficiently access rows and columns from a 10GB csv file which is too large to read into memory.
