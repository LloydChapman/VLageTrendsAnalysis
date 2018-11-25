clear
close all

%% Load processed data from running run_calc_CIs_and_ORs
load data

%% Set which studies to fit reversible catalytic model to
DATidx=find((strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (DAT)'))|strcmp(data.Author,'Koirala et al, 2004')|(strcmp(data.Author,'Ostyn et al, 2015')&strcmp(data.Type,'PREVALENCE'))|strcmp(data.Author,'Rijal et al, 2010')|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (DAT)'))|(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (DAT)')));
rK39idx=find(strcmp(data.Author,'Bern et al, 2007')|(strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (rK39)'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (rK39)')));
LSTidx=find(strcmp(data.Author,'Bern et al, 2006')|(strcmp(data.Author,'Nandy et al, 1987')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Patil et al, 2013')&strcmp(data.Type,'PREVALENCE (LST)'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (LST)'))|strcmp(data.Author,'Yangzom et al, 2012'));
PCRidx=find(strcmp(data.Author,'Kaushal et al, 2017')|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (PCR)')));

data.StartYear=str2double(cellfun(@(x) x(max(numel(x)-3,1):max(numel(x),0)),data.StartDate,'UniformOutput',false));

%% Fit const FOI model to each dataset
[p,p_test,p_all]=run_fit_Rev_Cat(data,DATidx,rK39idx,LSTidx,PCRidx,[],0,'RsltsAgeIndepFOI');

% Calculate summary statistics for estimated parameters for DAT studies
[med_lambda,med_gamma,age_pc_infctd,age_first_infctn,dur_imm_rspnse]=calc_smmry_stats('RsltsAgeIndepFOI.csv');

%% Fit const FOI model with same reversion rate across all studies
% DATidx=find((strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (DAT)'))|strcmp(data.Author,'Koirala et al, 2004')|strcmp(data.Author,'Rijal et al, 2010')|(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'PREVALENCE')));
[p_same_rvsn,p_same_rvsn_test,p_same_cnvsn_rvsn_test]=run_fit_Rev_Cat_same_rvsn(data,DATidx,rK39idx,LSTidx,PCRidx,[],0,'RsltsAgeIndepFOISameRvsn'); % N.B. Doesn't converge so no CIs if PCR data is included or Ostyn, Schenkel or Topno DAT data are included

%% Fit Hasker data
% Hasker seroprevalence and seroconversion data
dataH = [2 10 20 30 40 50 60 70 ;
        9 19 29 39 49 59 69 90 ;
        3858 2802 1565 1459 1021 812 767 321 ;
        100 126 92 123 101 87 97 51;
        3415 2146 1093 1084 765 602 555 213 ;
        65 59 54 44 45 38 31 15];
str = [];%'Hasker et al, 2013';
% Fit const FOI model to Hasker data
[par_c,parCI_c,NLL_c] = fit_Rev_Cat_Hasker(dataH,str,true,[],0);
% Fit age-dep FOI model to Hasker data
[par_a,parCI_a,NLL_a] = fit_Rev_Cat_Hasker(dataH,str,true,[],[]);

AIC_c=AIC(NLL_c,numel(par_c));
AIC_a=AIC(NLL_a,numel(par_a));
AICdiff=AIC_c-AIC_a;

%% Fit age-dep FOI model to each dataset with fitted slope from Hasker data
% paH=run_fit_Rev_Cat(data,DATidx,rK39idx,LSTidx,PCRidx,[],par_a(2),'RsltsAgeDepFOI_HaskerSlope.csv');

%% Fit age-dep FOI model to each dataset with intercept and slope (and gamma) as free parameters
[p_a,p_a_test]=run_fit_Rev_Cat(data,DATidx,rK39idx,LSTidx,PCRidx,[],[],'RsltsAgeDepFOI');

%% Fit const FOI model with same reversion rate across all studies/tests
[p_a_same_rvsn,p_a_test_same_rvsn,p_a_same_cnvsn_rvsn_test]=run_fit_Rev_Cat_same_rvsn(data,DATidx,rK39idx,LSTidx,PCRidx,[],[],'RsltsAgeDepFOISameRvsn');
