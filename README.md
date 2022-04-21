# Refactored Version of [scAge](https://github.com/alex-trapp/scAge)
Refactored version of Alex Trapp's [scAge](https://github.com/alex-trapp/scAge).

## Overview of Epigenetic Clocks

### CpG Sites
CpG sites are sites along the genome in which a cytosine is immediately followed by a guanine. In this pair, the cytosine can be methylated, meaning that a methyl group can be added to it. A CpG site can be in a methylated or unmethylated state, and there are processes which can methylate an unmethylated site and processes which can demethylate a methylated site. The state of a CpG site's methylation can affect the gene expression of the surrounding DNA sequence. The CpG sites themselves--meaning their existence and location in the genome--remain consistent across an organism's lifespan, and two organisms of the same species will have essentially identical CpG sites, just as they have essentially identical genomes. The methylation state of a CpG site can be different between two organisms, within the same organism at two different times, within the same organism in two different cells from two different tissues, or even within the same organism in two different cells in the same tissue. The overall methylation state of these CpG sites in an organism is called the methylome.

### Predicting Age
The methylome exhibits predictable changes with age. These age-related changes have been used to produce "epigenetic clocks". Trained on methylation data paired with actual "chronological" ages of samples, epigenetic clocks can take in methylation data from a handful of CpG sites and predict the age of a given sample. This age is interpreted as a "biological" age. In contrast to "chronological" age which represents the actual amount of time an organism has spent on this earth, "biological" age represents an organisms point along the underlying biolical process. If someone's biological age is below their chronological age, they are considered to have a lower rate of aging. If their biological age is above their chronological age, they are considered to have a higher rate of aging.

### Multi-Cell vs Single-Cell
