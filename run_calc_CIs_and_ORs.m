clear
close all
delete *.mat

%% Load database
data=readtable('S1_Data.xlsx');

%% Names of age group bounds variables
AgeLow={'Age1Low','Age2Low','Age3Low','Age4Low','Age5Low','Age6Low','Age7Low','Age8Low','Age9Low','Age10Low'};
AgeHigh={'Age1High','Age2High','Age3High','Age4High','Age5High','Age6High','Age7High','Age8High','Age9High','Age10High'};

%% Set which studies to plot
% DAT
% DATidx=find((strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (DAT)'))|strcmp(data.Author,'Koirala et al, 2004')|(strcmp(data.Author,'Ostyn et al, 2015')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Patil et al, 2013')&strcmp(data.Type,'PREVALENCE (DAT)'))|strcmp(data.Author,'Rijal et al, 2010')|(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (DAT)'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (DAT)')));
DATidx=find((strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (DAT)'))|strcmp(data.Author,'Koirala et al, 2004')|(strcmp(data.Author,'Ostyn et al, 2015')&strcmp(data.Type,'PREVALENCE'))|strcmp(data.Author,'Rijal et al, 2010')|(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (DAT)'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (DAT)')));
% rK39
rK39idx=find(strcmp(data.Author,'Bern et al, 2007')|(strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (rK39)'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (rK39)')));
% LST
LSTidx=find(strcmp(data.Author,'Bern et al, 2006')|(strcmp(data.Author,'Nandy et al, 1987')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Patil et al, 2013')&strcmp(data.Type,'PREVALENCE (LST)'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (LST)'))|strcmp(data.Author,'Yangzom et al, 2012'));
% PCR
PCRidx=find(strcmp(data.Author,'Kaushal et al, 2017')|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (PCR)')));
% Multiple tests
Brnidx=find((strcmp(data.Author,'Bern et al, 2007')&strcmp(data.Type,'PREVALENCE (rK39) (2002)'))|strcmp(data.Author,'Bern et al, 2006'));
Hskridx=find(strcmp(data.Author,'Hasker et al, 2013')&(strcmp(data.Type,'PREVALENCE (DAT)')|strcmp(data.Type,'PREVALENCE (rK39)')));
Sklidx=find(strcmp(data.Author,'Schenkel et al, 2006'));
Tpnidx=find(strcmp(data.Author,'Topno et al, 2010')&~strcmp(data.Type,'PREVALENCE'));
Mltidx={Brnidx,Hskridx,Sklidx,Tpnidx};

% VL
% Index for VL incidence studies with case numbers and study populations
VLidx=find(strcmp(data.Author,'Barnett et al, 2005')|(strcmp(data.Author,'Bern et al, 2005')&strcmp(data.Type,'INCIDENCE (RAW DATA)'))|strcmp(data.Author,'Ferdousi et al, 2012')|strcmp(data.Author,'Hasker et al, 2012')|(strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'INCIDENCE'))|strcmp(data.Author,'Picado et al, 2014')|(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'INCIDENCE')));

%% Calculate confidence intervals (CIs) and odds ratios (ORs)
[data,SeroTests,SeroPoss,SeroPosPrevs,CIL,CIU,OR,ORCIL,ORCIU,p,StudyPops,VLCases,VLIncs,VLCIL,VLCIU,VLOR,VLORCIL,VLORCIU,VLp,VLRR,VLRRCIL,VLRRCIU,VLRRp,VLIRR,VLIRRCIL,VLIRRCIU,VLIRRp]=calc_CIs_and_ORs(data);

%% Run chi-square trend test on VL incidence for age>20
data=run_chi_sq_trend_test(data,VLidx,AgeLow,AgeHigh,StudyPops,VLCases);

%% Extract raw data for each study from database
extract_raw_data(data,AgeLow,AgeHigh,SeroTests,SeroPoss)

%% Plot seroprevalence and VL incidence by age with CIs
plot_seroprev_and_VL_inc_with_CIs(data,AgeLow,AgeHigh,DATidx,rK39idx,LSTidx,PCRidx,Mltidx,SeroPosPrevs,CIL,CIU,VLidx,VLIncs,VLCIL,VLCIU)

%% Print out table for each study
print_table(data,AgeLow,AgeHigh,SeroTests,SeroPoss,SeroPosPrevs,CIL,CIU,OR,ORCIL,ORCIU,p,StudyPops,VLCases,VLIncs,VLCIL,VLCIU,VLRR,VLRRCIL,VLRRCIU,VLRRp)


