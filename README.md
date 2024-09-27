# Gadea-Gonzalo-Ramos-Trends-2024

This README file clarifies the MATLAB codes to replicate the paper "Trends in Temperature Data: Micro-foundations of Their Nature", by Gadea, M.D., Gonzalo, J. and A. Ramos.

The folders EMPIRICAL and SIMULATIONS contain the codes to reproduce the results of the paper.

The folder FUNCTIONS contains several auxiliary funcions

REMEMBER CHANGING THE PATH ACCORDINGLY

-------------------------------------------------------------------------------------------------
EMPIRICAL
The code 'read_grid_file.m' reads the file 'CRUTEM.5.0.1.0.anomalies.nc' which can be downloaded from 
https://www.metoffice.gov.uk/hadobs/crutem5/data/CRUTEM.5.0.1.0/download.html and generates de data 
structure 'DATA_GRIDS.mat'

The code 'make_figure1.m' produces Figures 1a and 1b of the paper.
The code 'calculate_grids_selected.m' computes the selected grids and buil the data structure 
'grids_yearly_selected_1880_2022.mat' or 'grids_yearly_selected_1960_2022.mat'. You must properly change 
the options in the file to get both outputs.

The code 'main_grids_1880_2022.m', construct the global, NH, SH temperature averages with method A and B.
It generates Figure 2 of the paper and computes the ADF test that are displayed in Table 1. The average are
save in the data structure 'DATA_GRID_1880_2022.mat'.

The code 'main_UR_SB_1880_2022.m', computes the Kim-Perron-2009 tests for unit root with structural breaks;
the results are displayed in Table 1.

The code 'main_structural_breaks_1880_2022.m' computes the Perron-Yabu-2009 test of structural breaks for mean
method A, method B and individual series; the results are displayed in Table 1.

The code 'main_Table2.m' reproduces the Table 2 of the paper and produces Figures 3a and 3b.



SIMULATIONS
The code 'main_grids_compute_parameters.m' reads the data structure 'DATA_GRID.mat' and:
0.- Test for unit root
1.- Estimate a linear trend model
2.- Estimate a cuadratic trend model
3.- Estimate a cubic trend model
4.- Detect the optimal position of a break in constant and trend and estimate the model
5.- Estimate a linear trend model since 1960, with observations present in the complete sample, 1880-2022, since 1960
6.- Only observations from 1960; we work with observations present in the period 1960-2022
7.- With observations since 1960; we work with observations that are new in the period 1960-2022

and generates the structure 'PARAM_GRIDS.mat'

The main script is 'simulations_EL' that carries out the simulations, saves the results in the structure RDOS_simus
The scripts 'analysis_rdos.m' and 'build_table_latex.m' analyses the results and produces the Latex Tables 3 and 4 of the paper



Good luck,
Lola, Jesus and Andrey

PS: Remember, any code works at the first try.
