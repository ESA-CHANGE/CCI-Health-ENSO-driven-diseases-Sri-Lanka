# Case of Study Title

ENSO-Driven Climate Variability and Dengue & Leptospirosis in Sri Lanka

## Description

This case study investigates how large-scale climate variability associated with the El Niño–Southern Oscillation (ENSO) influences dengue and leptospirosis risk across Sri Lanka, and how Earth observation–derived environmental conditions represented by ESA Climate Change Initiative (ESA CCI) Essential Climate Variables (ECVs) contribute to spatial differences in climate-related disease risk. District-level monthly disease surveillance data are analysed alongside ENSO indicators, including the Oceanic Niño Index (ONI) and Dipole Mode Index (DMI), and environmental observations such as rainfall and vegetation dynamics derived from satellite-based Earth observation data.
ESA CCI ECV products provide standardized, spatially consistent environmental information on land cover, vegetation, rainfall, soil moisture, and human settlement characteristics, which are used to characterize the environmental context of each district. These ECV-derived indicators enable assessment of whether differences in land-surface and environmental conditions modify or explain spatial heterogeneity in ENSO–disease associations.
The epidemiological analysis uses a two-stage modelling framework. District-specific distributed-lag models first estimate delayed associations between ENSO variability and disease incidence. These estimates are subsequently combined using multivariate random-effects meta-analysis, followed by meta-regression incorporating ESA CCI ECV-derived environmental and land-surface characteristics to identify potential drivers of spatial heterogeneity.
The overall objective is to demonstrate how ESA ECVs can be integrated with epidemiological surveillance and climate indicators to quantify climate-related health risks, characterize spatial vulnerability, and support climate-health surveillance, risk assessment, and adaptation planning in Sri Lanka.
Data
Epidemiological and Climate Time-Series Data

The analysis_df dataset contains district-month observations for Sri Lanka and forms the main longitudinal dataset for the epidemiological analysis. It includes:

District and province identifiers
District population and P-codes
Year and month
Monthly dengue cases
Monthly leptospirosis cases
Oceanic Niño Index (ONI)
Dipole Mode Index (DMI)
Monthly rainfall totals
Rainfall relative to the long-term average
Monthly and rolling 1-month and 3-month rainfall measures
Rainfall anomaly indicators
Monthly mean NDVI
NDVI long-term averages and anomalies
Niño 3.4 sea-surface temperature and anomaly measures

These variables provide the temporal exposure and outcome series used in the first-stage distributed-lag analyses.

District-Level Environmental Covariates

The covariate_matrix dataset contains district-level environmental and Earth observation-derived characteristics used in the second-stage analysis. Variables include:
Land-cover composition based on LCCS classes
Dominant land-cover class
Human Settlement Index (HSI) for 2016 and 2025
Absolute and relative changes in HSI
Long-term mean rainfall
Annual district-level rainfall characteristics
Long-term mean NDVI
Annual NDVI characteristics
Changes in rainfall and NDVI between 2007 and 2025
Annual soil-moisture estimates
Long-term mean soil moisture
Absolute and relative changes in soil moisture
Number of valid years contributing to soil-moisture estimates

These variables are used to investigate potential explanations for between-district heterogeneity in the climate–disease associations.

## Getting Started

### Dependencies

R: version >= 4.3
Operating system: Windows, macOS, or Linux
Required R packages:
dplyr
dlnm
mvmeta
splines
ggplot2

dplyr is used for data preparation and manipulation, dlnm for distributed-lag modelling, mvmeta for multivariate random-effects meta-analysis and meta-regression, splines for flexible exposure-response and temporal functions, and ggplot2 for visualization.

Additional packages may be required by specific data-preparation or Earth observation processing scripts.

### Installing

Install the required R packages using:

install.packages(c(
  "dplyr",
  "dlnm",
  "mvmeta",
  "splines",
  "ggplot2"
))

The splines package is distributed with R and normally does not require separate installation.

### Executing program

The analysis follows a sequential two-stage workflow.

Step 1 – First-stage distributed-lag models

District-specific distributed-lag models are fitted to estimate the association between ONI and disease incidence over biologically plausible lag periods.

source("code/Case_6_Model_Lepto_1.R")

The first-stage models generate district-specific coefficient estimates and their corresponding variance-covariance matrices.

Step 2 – Second-stage meta-analysis

The district-specific estimates are combined using a multivariate random-effects meta-analysis to obtain pooled associations and quantify between-district heterogeneity.

source("code/Case_6_Model_Lepto_1.R")

Step 3 – Univariable meta-regression

District-level environmental and Earth observation covariates are evaluated individually to investigate potential explanations for heterogeneity in the ONI–disease association.

source("code/Case_6_Model_Lepto_2.R")

Covariates include land-cover characteristics, human settlement, rainfall, NDVI, and soil moisture.

Step 4 – Multivariable meta-regression

A scientifically selected set of district-level environmental covariates is included simultaneously in a multivariable meta-regression to assess whether combinations of environmental characteristics explain between-district heterogeneity.

source("code/Case_6_Model_Lepto_2.R")

Step 6 – Sensitivity analyses

Sensitivity analyses evaluate the robustness of the estimated climate–disease associations to alternative lag specifications, model assumptions, and environmental covariate definitions.

source("code/Case_6_Model_Lepto_Sensitivity.R")

```
code blocks for commands
```

## Help

Any advise for common problems or issues.

```
command to run if program contains helper info
```

## Authors
Prasad Liyanage (PhD)


## Version History

Version 1.0 – Initial release of the ENSO–Leptospirosi analytical workflow.
Development version – Data preparation, distributed-lag modelling, multivariate meta-analysis, meta-regression, and sensitivity analyses.

## License

This project is licensed under the MIT License of 2023 ESA Climate Change Initiative - see the LICENSE file for details.

## Acknowledgments

This work was developed within the Climate-Health Adaptation Through New Generation Earth Observations (CHANGE) project and acknowledges the contribution of the European Space Agency and collaborating institutions.
We acknowledge the Sri Lankan health authorities and surveillance systems contributing to the underlying epidemiological data, as well as the providers of the climate, Earth observation, land-cover, vegetation, rainfall, and soil-moisture datasets used in this analysis.
