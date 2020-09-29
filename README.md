# Taxon-Specific Phytoplankton Growth, Nutrient Utilization, and Light Limitation in the Oligotrophic Gulf of Mexico
*NATALIA YINGLING, THOMAS B. KELLY, TAYLOR A. SHROPSHIRE, MICHAEL R. LANDRY, KAREN E. SELPH, ANGELA N. KNAPP, SVEN KRANZ, MICHAEL R. STUKEL*
### About the Project

Using data collected on two field expeditions into the oligotrophic Gulf of Mexico (GoM), we use an objective optimization procedure to determine a *NEMURO-like* parameterization for 5 distinct phytoplankton taxa (PRO, SYN, DIA, ADINO, and PRYM) and one residual group (i.e. OTHER). Doing so provided both an optimum parameter set as well as a collection of distinct parameterizations with similar objective cost (i.e. log-likelihood), which could be tested in a full, 3D biogeochemical domain.

Data-assimilation was conducted through a Markov Chain Monte Carlo approach whereby the original NEMURO-GoM parameter set (Shropshire et al. 2019) was randomly perturbed and conditionally accepted based on model-observation misfit and prior likelihood. This procedure iterated > 10<sup>6</sup> times until a stationary distribution was reached. The resulting solution sets provided a distribution for each parameter (e.g. half-saturation coefficient).

All scripts and data necessary to run this model are provided as part of this REPO.

#### Running the model
1. First, a recent version of R is required. We recommend RStudio as a convenient front-end and IDE.
2. Download the files in this REPO: Green code button > Download Zip.
3. Included are the RStudio .Rproj file, this readme, and three directories:
 *  Data: The data directory contains the xlsx spreadsheets and raw data used in the model.
 * R: The R folder contains all r scripts to run and analyze model output.
 * Yingling et al. Manuscript: Folder contains a preprint and figures of the model as well as the final dataset (.rdata format) and sensitivity tests.
4. The R scripts (inside the R folder) are split into two files: *MCMC Functions.R* which contains the functions used to run the model, and *MCMC Project.R* which calls the functions and is the user interactive side of the model.
5. Simply run **MCMC Project.R** to recreate the dataeset presented (including model figures).


<hr />

### Abstract & Citation

__ABSTRACT__
The vast oligotrophic regions of the oceans are predominantly nitrogen limited in the surface ocean and light limited at the deep chlorophyll maximum (**DCM**); hence determining light and nitrogen co-limitation patterns for diverse phytoplankton taxa is crucial to understanding marine primary production. On 2 month-long cruises (May 2017 and 2018) in the open-ocean Gulf of Mexico, we measured primary productivity, nitrate uptake, and ammonium uptake throughout the euphotic zone. Primary productivity generally declined with depth from the mixed layer to the DCM, with vertically-integrated values that averaged 27.1 mmol C m<sup>-2</sup> d<sup>-1</sup>.  *f*-ratios (=nitrate uptake/nitrate uptake+ammonium uptake) were consistently low with average upper euphotic zone values ranging from 0.01 – 0.14 and lower euphotic zone values ranging from 0.03 – 0.44.  We used a Markov Chain Monte Carlo statistical approach to assimilate rate measurements and taxon-specific phytoplankton biomasses assessed at the same locations and depth into a model of taxon-specific phytoplankton nutrient and light-limitation patterns. We found that the field data here indicates that the oligotrophic GoM is a region that is pico-plankton dominated, mostly by Prochlorococcus, highly stratified with low subsurface nitrate and chl a. concentration with a deep DCM, typically 80-115 m. The data suggests that this region strongly depends on recycled production due to the fact ammonium uptake is on the scale of 2-3 fold higher than nitrate uptake.

__Citation__

Yingling, N. et al. Taxon-Specific Phytoplankton Growth, Nutrient Limitation, and Light Limitation in The Oligotrophic Gulf of Mexico. Journal of Plankton Research.
