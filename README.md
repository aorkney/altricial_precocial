This repository contains R scripts and data required to replicate
the analyses in Orkney et al. 2025. 


Growing up and fledgling is an exhausting challenge for chicks, which may restrict the variety of adult forms birds can evolve.  Many bird species assist their highly dependent chicks, providing food, warmth and shelter, and we hypothesised that this might alleviate developmental restrictions in these lineages in contrast to species with less supportive parenting. We compared adult limb bone length and flight-style descriptions across a representative sample of species. We fit evolutionary models, showing that parental support is associated with faster evolutionary divergences between the wing and leg skeleton and the exploration of a wider variety of flight-style combinations.


In all analyses we employ the all-birds super tree recently published by 
McTavish et al., 2025. (https://doi.org/10.1073/pnas.2409658122)

Data:

Ontogeny_Bjarnason_04_30_2025.csv
This spreadsheet contains limb bone length measurements derived from 
osteological measurements across a variety of bird species from
Bjarnason et al. 2021 (https://doi.org/10.18563/journal.m3.125), 
body mass estimates assembled from Cornell Birds of the World (https://birdsoftheworld.org/bow/support/citations-and-references),
 and categorical designations for developmental strategy sourced from 
Starck et al. 1993. (https://doi.org/10.1007/978-1-4615-9582-3_6)

ontogeny_magnificent.csv
This spreadsheet contains a small subset limb bone length measurements aggregated by 
Brinkworth et al. 2023 (https://doi.org/10.1038/s41467-023-41415-2), and associated metadata including body mass estimates.
I have appended categorical designations for developmental strategy sourced from 
Starck et al. 1993. (https://doi.org/10.1007/978-1-4615-9582-3_6)

flight_masses_22_10_2022_plus_A.csv
This spreadsheet contains binary flight-style designations for the species in 
Bjarnason et al. 2021 (https://doi.org/10.18563/journal.m3.125), 
compiled as an extension of (Taylor, G.K. and Thomas, A.L., 2014. Evolutionary biomechanics: selection, phylogeny, and constraint. Oxford University Press (UK) ),
and originally published in Orkney et al. 2021
(https://doi.org/10.1038/s41559-021-01509-w) 

Scripts:

Figure_1_04_08_2025.R 
This script produces the analysis for Figure 1. 
Body mass dependent scaling of limb bone lengths is removed from limb bone length measurements,
following a linear regression approach in which non-independence between species
is represented by their shared ancestry and a random walk model of evolution is assumed.
Thereafter the Pearson's correlations are computed between pairs of the residuals in cohorts of
bird species with Altricial and Preoccial developmental strategies. 
We apply a Fisher-transform and test the one-tailed hypothesis that correlation magnitudes
are lower in Altricial cohorts, accommodating sample non-independence under the assumption
of a random walk model of evolution. 

We show that the wing and leg proportions of Altricial birds are significantly less correlated
over their evolutionary history than are those of Precocial birds. 
This conforms to our hypothesis that the high parental investment of Altricial developmental
strategies alleviates physiological and developmental restrictions that otherwise
reduce the capacity for the independent evolution of wing and leg proportions in birds. 

See the in-line comments for specific technical descriptions and citations. 

Figure_2_04_08_2025.R 
This script produces a phenogram/'trait gram' that describes how Altricial and Precocial
cohorts of birds have explored divergent evolution of their wing and leg proportions over their evolutionary history. Moreover the relationship that divergent wing and leg evolution has
to ecological adaptation and body mass across Altricial and Precocial cohorts is qualitatively explored. 
We first take residual limb proportions and compute a singular value decomposition on the covariance matrix between wing and leg residuals. 
The first axis explains the vast majority of variance, and we therefore treat this as a vector
of concordant evolution between the wing and leg. 
The orthogonal complement to this vector is then taken as an axis of divergence between 
wing and leg. 
We compute the divergence of all bird species by projecting their residuals upon this vector. 
Thereafter, we compute a stochastic character map of Altricial and Precocial developmental
strategies across birds, assuming Precocial at the root and equal transition rates. 
Thereafter, we compute ancestral state reconstructions of wing-leg divergence across birds, using the stochastic character map to partition our analysis to resolve different Brownian rates of variance. We conduct a likelihood ratio test comparing this model to a single-rate
model. 
We also take avian flight-style categorisations, perform a Gower transform and find the 
principal coordinate axes. We then find the Mahalanobis distances of individual birds
from the mean position in these axes. 
We fit polynomial regressions in Altricial and Precocial cohorts of bird species depending
on wing-leg divergence to investigate whether divergent wing and leg evolution may afford
access to more unusual flight-style ecologies and whether developmental strategy might
moderate this process. 
We also fit polynomial regressions in Altricial and Precocial cohorts of bird species
where body mass is the dependent variable, allowing us to explore whether ontogeny/growth trajectory facilitates wing-leg divergence and if this differs depending on avian 
developmental strategy. 

We show that Altricial birds explore a wider variety of wing-leg divergence than do
Precocial bird lineages, and that ancestral state reconstructions show this pattern is 
consistent with a faster rate of evolutionary divergence between these structures in Altricial birds. 
Moreover, we show that highly divergence wing and leg proportions in Altricial birds are associated with unusual flight-style combinations and that variety in wing-leg divergence in Altricial birds is significantly associated with body mass variety. This implies
that a growth-dependent mechanism facilitated wing-leg evolutionary divergence in Altricials
and that this process is materially relevant to their ecological adaptation. 

See the in-line comments for specific technical descriptions, package dependences and citations. 

Figure_3_04_08_2025.R 
This figure produces illustrative ternary diagrams of internal wing proportions, 
colourised by discretised developmental categorisations from Starck et al., 1993 
(https://doi.org/10.1007/978-1-4615-9582-3_6), and by the estimated divergence between
wing and leg evolution.
These reveal that there are two major axes of variation in internal wing proportions in
birds. Most birds, and especially those with relatively low divergence between their wing and
leg skeletal traits, conform to the primary axis of variation. 
Variation perpendicular to this axis (associated with 'brachial index', a kinematically 
relevant wing characteristic in birds) is associated with more divergent wing and leg 
skeletal traits and with developmental categorisations that are more Altricial in character. 


