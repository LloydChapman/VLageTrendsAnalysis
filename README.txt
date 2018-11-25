Title: MATLAB code for performing statistical analyses and fitting reversible catalytic models in 'Age trends in asymptomatic and symptomatic Leishmania donovani infection in the Indian subcontinent: a review and analysis of data from diagnostic and epidemiological studies'

Version: 0.1.0

Author: Lloyd Chapman

Email: Lloyd.Chapman@lshtm.ac.uk

Description: MATLAB code for performing statistical analyses and fitting age-independent and age-dependent reversible catalytic models described in [1], using age-stratified data on asymptomatic Leishmania donovani infection prevalence and clinical VL incidence from different studies in the Indian subcontinent (S1_Data). Further details on the models and model fitting approach are provided in [2]. The main files for performing the statistical analyses and fitting the catalytic models are "run_calc_CIs_and_ORs.m" and "gen_Rev_Cat_rslts.m" respectively. Running these files sequentially produces Figs 2-5 in [1] and Figs 1-3 and Tables 2-4 in [2]. The code uses the MATLAB File Exchange files "cochran_arm.m" and "distinguishable_colors.m", which are available from https://uk.mathworks.com/matlabcentral/fileexchange/45186-cochran-armitage-test and https://uk.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors, and "confint.m", which is available at https://brainder.org/2012/04/21/confidence-intervals-for-bernoulli-trials/.

Repository contents:
add_total.m
age_dep_FOI.m
AIC.m
calc_CIs_and_ORs.m
calc_smmry_stats.m
cochran_arm.m
cochran_arm_license.txt
confint.m
distinguishable_colors.m
extract_raw_data.m
fit_Rev_Cat.m
fit_Rev_Cat_Hasker.m
fit_Rev_Cat_same_rvsn.m
gen_Rev_Cat_rslts.m
negLL_Hasker.m
negLL_ODE.m
negLL_ODE0.m
negLL_ODE_same_rvsn.m
plot_cnvsn_vs_rvsn_rates.m
plot_LL.m
plot_LL_ODE.m
plot_model_fit.m
plot_seroprev_and_VL_inc_with_CIs.m
poiss_CI.m
print_par_ests.m
print_table.m
put_in_table.m
README.txt
Rev_Cat_ODE.m
run_calc_CIs_and_ORs.m
run_chi_sq_trend_test.m
run_fit_Rev_Cat.m
run_fit_Rev_Cat_same_rvsn.m
solve_Rev_Cat_ODE.m
S1_Data.xlsx
(See individual files for a description of their function.)

Developed in: MATLAB R2016b (9.1.0) @ 1984-2016 The MathWorks, Inc. 

Installation: MATLAB, which requires a user license, must be installed to run the code. It can be downloaded from https://uk.mathworks.com/downloads/. Information on installing and activating MATLAB can be found at https://uk.mathworks.com/help/install/. Save the contents of this archive in a folder, open MATLAB and change the current directory to the folder in which the files are saved. The code for plotting the data and calculating the odds ratios and risk ratios can then be run by typing “run_calc_CIs_and_ORs” at the MATLAB command prompt. The code for fitting the reversible catalytic models can be run by typing "gen_Rev_Cat_rslts" at the prompt (N.B. run_calc_CIs_and_ORs.m must be run first).

License: GNU Affero General Public License v3.0 (http://www.gnu.org/licenses/agpl-3.0.txt)

References: [1] Chapman LAC, et al. Age trends in asymptomatic and symptomatic Leishmania donovani infection in the Indian subcontinent: a review and analysis of data from diagnostic and epidemiological studies. PLoS Neglected Tropical Diseases. 2018; 12(12): e0006803. https://doi.org/10.1371/journal.pntd.0006803
[2] Chapman LAC, et al. Age trends in asymptomatic and symptomatic Leishmania donovani infection in the Indian subcontinent: a review and analysis of data from diagnostic and epidemiological studies - S2 Text. PLoS Neglected Tropical Diseases. 2018; 12(12): e0006803. https://doi.org/10.1371/journal.pntd.0006803
