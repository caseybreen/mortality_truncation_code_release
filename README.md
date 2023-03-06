## Mortality Modeling of Partially Observed Cohorts Using Administrative Death Records

This repository contains code and materials to replicate ["Mortality Modeling of Partially Observed Cohorts Using Administrative Death Records."](https://https://doi.org/10.31235/osf.io/efdzh)

### Replication Package

The this repository includes code to replicate all figures and tables in the paper. There are three steps to running the replication code: 

1. Clone this repository
2. Download the `data.zip` file from the [accompanying OSF project](https://doi.org/10.17605/OSF.IO/D6QHF), unzip the data file, and move it to the top level of the repository. 
3. Run the `00_run_all.Rmd` script, which will run all code (or run all scripts individually in any)


#### Data 

Please download all data for replication from the project's [OSF project](https://doi.org/10.17605/OSF.IO/D6QHF). The data were originally obtained from: 

- IPUMS-USA [[link](https://usa.ipums.org/usa/)]
- CenSoc [[link](https://censoc.berkeley.edu/)]

#### Code 

After downloading the required data and moving into the top level of the replication repository, researchers can run the following script to replicate all figures and tables: 

- `00_run_all.Rmd` - this file runs all scripts. 

Alternatively, researchers can run each script individually in any order. 

- `fig1_simulation_ols_window_width.Rmd` - generates Figure 1, which illustrates the bias of OLS regression under double truncation using a simulation study. 
- `fig2_regression_example.Rmd` - generates Figure 2, which provides an illustration of why regression on age of death can be biased under double truncation. 
- `fig3_sweden_example.Rmd` - generates Figure 3, which gives an example of estimated mortality differentials between men and women in Sweden (using HMD data) under truncation and not under truncation. 
- `fig4_censoc_dmf_1910_diagnostics.Rmd` - generates Figure 4, which applies graphical diagnostics from the Gompertztrunc package to a case study of the relationship between education and longevity. 
- `fig5_association_between_education_longevity.Rmd` - generates Figure 5, which shows the estimated association between education and longevity estimated from OLS regression and from our Gompertz MLE approach. 
- `fig6_simulation_window_numident.Rmd` - generates Figure 6, which simulates different windows widths of mortality observation and its effect on estimates from our Gompertz MLE approach. 

### Authors

- [Joshua R. Goldstein](https://jrgoldstein.com/)
- [Casey F. Breen](caseybreen.com)
- [Serge Atherwood](https://satherwood.wordpress.com/)
- Maria Osborne

