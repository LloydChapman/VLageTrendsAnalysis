function [data,SeroTests,SeroPoss,SeroPosPrevs,CIL,CIU,OR,ORCIL,ORCIU,p,StudyPops,VLCases,VLIncs,VLCIL,VLCIU,VLOR,VLORCIL,VLORCIU,VLp,VLRR,VLRRCIL,VLRRCIU,VLRRp,VLIRR,VLIRRCIL,VLIRRCIU,VLIRRp]=calc_CIs_and_ORs(data)
%CALC_CIS_AND_ORS Calculate confidence intervals and odds/risk ratios for
%prevalences and incidences

%% CALCULATE BINOMIAL CIs FOR SEROPREVALENCE
% Index for studies with seroprevalence and number tested but not numbers of seropositives
noSeroPos=find((strcmp(data.Author,'Hasker et al, 2013')&(strcmp(data.Type,'PREVALENCE (DAT)')|strcmp(data.Type,'PREVALENCE (rK39)')))|strcmp(data.Author,'Patil et al, 2013')|strcmp(data.Author,'Rijal et al, 2010'));
% Variables for calculating CIs
SeroTests={'SeroTest1','SeroTest2','SeroTest3','SeroTest4','SeroTest5','SeroTest6','SeroTest7','SeroTest8','SeroTest9','SeroTest10'};
SeroPoss={'SeroPos1','SeroPos2','SeroPos3','SeroPos4','SeroPos5','SeroPos6','SeroPos7','SeroPos8','SeroPos9','SeroPos10'};
SeroPosPrevs={'SeroPosPrev1','SeroPosPrev2','SeroPosPrev3','SeroPosPrev4','SeroPosPrev5','SeroPosPrev6','SeroPosPrev7','SeroPosPrev8','SeroPosPrev9','SeroPosPrev10'};
CIL={'CIL1','CIL2','CIL3','CIL4','CIL5','CIL6','CIL7','CIL8','CIL9','CIL10'};
CIU={'CIU1','CIU2','CIU3','CIU4','CIU5','CIU6','CIU7','CIU8','CIU9','CIU10'};
% Calculate binomial confidence intervals for prevalences
% Significance level
alpha=0.05;
zalpha=norminv(1-alpha/2);
% Method for calculating binomial confidence intervals
meth='Clopper-Pearson';
for i=1:numel(SeroTests)
    % Calculate missing seropositive numbers
    data.(SeroPoss{i})(noSeroPos)=round(data.(SeroTests{i})(noSeroPos).*data.(SeroPosPrevs{i})(noSeroPos));
    [data.(CIL{i}),data.(CIU{i})]=confint(data.(SeroTests{i}),data.(SeroPoss{i}),alpha,meth);
end
% Calculate missing total number seropositive
noSeroPosTot=find((strcmp(data.Author,'Hasker et al, 2013')&(strcmp(data.Type,'PREVALENCE (DAT)')|strcmp(data.Type,'PREVALENCE (rK39)')))|strcmp(data.Author,'Patil et al, 2013'));
data.SeroPos(noSeroPosTot)=sum(data{noSeroPosTot,SeroPoss},2,'omitnan');
data.SeroPosPrev(noSeroPosTot(1:2))=data.SeroPos(noSeroPosTot(1:2))./data.SeroTest(noSeroPosTot(1:2));

%% CALCULATE POISSON CIs FOR VL INCIDENCE
% Index for VL incidence studies with case numbers and study populations
StudyPops={'StudyPop1','StudyPop2','StudyPop3','StudyPop4','StudyPop5','StudyPop6','StudyPop7','StudyPop8','StudyPop9','StudyPop10'};
VLCases={'VLCase1','VLCase2','VLCase3','VLCase4','VLCase5','VLCase6','VLCase7','VLCase8','VLCase9','VLCase10'};
VLIncs={'VLInc1','VLInc2','VLInc3','VLInc4','VLInc5','VLInc6','VLInc7','VLInc8','VLInc9','VLInc10'};
VLCIL={'VLCIL1','VLCIL2','VLCIL3','VLCIL4','VLCIL5','VLCIL6','VLCIL7','VLCIL8','VLCIL9','VLCIL10'};
VLCIU={'VLCIU1','VLCIU2','VLCIU3','VLCIU4','VLCIU5','VLCIU6','VLCIU7','VLCIU8','VLCIU9','VLCIU10'};
% Index for studies without case numbers
noVLCase=strcmp(data.Author,'Barnett et al, 2005');
% Calculate Poisson confidence intervals for incidence rates
for i=1:numel(VLCases)
    % Calculate missing case numbers
    data.(VLCases{i})(noVLCase)=round(data.(VLIncs{i})(noVLCase).*data.(StudyPops{i})(noVLCase).*data.NumYrVLInc(noVLCase)/1000);
    % Recalculate incidences from calculated case numbers
    data.(VLIncs{i})(noVLCase)=data.(VLCases{i})(noVLCase)/(data.(StudyPops{i})(noVLCase)*data.NumYrVLInc(noVLCase))*1000;
    [data.(VLCIL{i}),data.(VLCIU{i})]=poiss_CI(data.(VLCases{i}),data.(StudyPops{i}).*data.NumYrVLInc/1000,alpha);
end
% Enter CIs for Singh et al 2010 manually
noVLCaseOrStudyPop=find(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'INCIDENCE'));
VLCILSingh=[0.89,6.57,4.32,4.37,2.56,0.96,NaN(1,4)];
VLCIUSingh=[3.90,10.52,7.63,8.38,7.56,5.51,NaN(1,4)];
for i=1:numel(VLCIL)
data.(VLCIL{i})(noVLCaseOrStudyPop)=VLCILSingh(i);
data.(VLCIU{i})(noVLCaseOrStudyPop)=VLCIUSingh(i);
end

%% CALCULATE ODDS RATIOS FOR SEROPREVALENCE WITH AGE
SeroNeg={'SeroNeg1','SeroNeg2','SeroNeg3','SeroNeg4','SeroNeg5','SeroNeg6','SeroNeg7','SeroNeg8','SeroNeg9','SeroNeg10'};
OR={'OR1','OR2','OR3','OR4','OR5','OR6','OR7','OR8','OR9','OR10'};
ORCIL={'ORCIL1','ORCIL2','ORCIL3','ORCIL4','ORCIL5','ORCIL6','ORCIL7','ORCIL8','ORCIL9','ORCIL10'};
ORCIU={'ORCIU1','ORCIU2','ORCIU3','ORCIU4','ORCIU5','ORCIU6','ORCIU7','ORCIU8','ORCIU9','ORCIU10'};
p={'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10'};

% Re-order data if no seropositives in youngest age group
if sum(data.SeroPos1==0)~=0
    zidx=(data.SeroPos1==0);
    l=find(table2array(data(zidx,SeroPoss))~=0,1);
    m=find(~isnan(table2array(data(zidx,SeroPoss))),1,'last');
    nzidx=[l:m,1:l-1,m+1:numel(SeroTests)];
    data(zidx,SeroTests)=data(zidx,SeroTests(nzidx));
    data(zidx,SeroPoss)=data(zidx,SeroPoss(nzidx));
    data(zidx,SeroPosPrevs)=data(zidx,SeroPosPrevs(nzidx));
    reord=true;
end

data.SeroNeg1=data.SeroTest1-data.SeroPos1;
n=size(data,1);
data.OR1=NaN(n,1);
data.ORCIL1=NaN(n,1);
data.ORCIU1=NaN(n,1);
data.p1=NaN(n,1);
for i=2:numel(SeroTests)
    data.(SeroNeg{i})=data.(SeroTests{i})-data.(SeroPoss{i});
    data.(OR{i})=(data.(SeroPoss{i}).*data.SeroNeg1)./(data.(SeroNeg{i}).*data.SeroPos1);
    stderr=sqrt(1./data.(SeroPoss{i})+1./data.(SeroNeg{i})+1./data.SeroPos1+1./data.SeroNeg1);
    z=log(data.(OR{i}))./stderr;
    data.(p{i})=2*normcdf(-abs(z));
    data.(ORCIL{i})=exp(log(data.(OR{i}))-zalpha*stderr);
    data.(ORCIU{i})=exp(log(data.(OR{i}))+zalpha*stderr);
end

% Re-order data so numbers (OR etc.) match age groups
if reord
    rnzidx=[(m-l+2):m,1:(m-l+1),(m+1):10];
    data(zidx,SeroTests)=data(zidx,SeroTests(rnzidx));
    data(zidx,SeroPoss)=data(zidx,SeroPoss(rnzidx));
    data(zidx,SeroPosPrevs)=data(zidx,SeroPosPrevs(rnzidx));
    data(zidx,SeroNeg)=data(zidx,SeroNeg(rnzidx));
    data(zidx,OR)=data(zidx,OR(rnzidx));
    data(zidx,p)=data(zidx,p(rnzidx));
    data(zidx,ORCIL)=data(zidx,ORCIL(rnzidx));
    data(zidx,ORCIU)=data(zidx,ORCIU(rnzidx));
end

%% CALCULATE ODDS RATIOS FOR VL INCIDENCE WITH AGE
NonVLCases={'NonVLCase1','NonVLCase2','NonVLCase3','NonVLCase4','NonVLCase5','NonVLCase6','NonVLCase7','NonVLCase8','NonVLCase9','NonVLCase10'};
VLOR={'VLOR1','VLOR2','VLOR3','VLOR4','VLOR5','VLOR6','VLOR7','VLOR8','VLOR9','VLOR10'};
VLORCIL={'VLORCIL1','VLORCIL2','VLORCIL3','VLORCIL4','VLORCIL5','VLORCIL6','VLORCIL7','VLORCIL8','VLORCIL9','VLORCIL10'};
VLORCIU={'VLORCIU1','VLORCIU2','VLORCIU3','VLORCIU4','VLORCIU5','VLORCIU6','VLORCIU7','VLORCIU8','VLORCIU9','VLORCIU10'};
VLp={'VLp1','VLp2','VLp3','VLp4','VLp5','VLp6','VLp7','VLp8','VLp9','VLp10'};
VLRR={'VLRR1','VLRR2','VLRR3','VLRR4','VLRR5','VLRR6','VLRR7','VLRR8','VLRR9','VLRR10'};
VLRRCIL={'VLRRCIL1','VLRRCIL2','VLRRCIL3','VLRRCIL4','VLRRCIL5','VLRRCIL6','VLRRCIL7','VLRRCIL8','VLRRCIL9','VLRRCIL10'};
VLRRCIU={'VLRRCIU1','VLRRCIU2','VLRRCIU3','VLRRCIU4','VLRRCIU5','VLRRCIU6','VLRRCIU7','VLRRCIU8','VLRRCIU9','VLRRCIU10'};
VLRRp={'VLRRp1','VLRRp2','VLRRp3','VLRRp4','VLRRp5','VLRRp6','VLRRp7','VLRRp8','VLRRp9','VLRRp10'};
VLIRR={'VLIRR1','VLIRR2','VLIRR3','VLIRR4','VLIRR5','VLIRR6','VLIRR7','VLIRR8','VLIRR9','VLIRR10'};
VLIRRCIL={'VLIRRCIL1','VLIRRCIL2','VLIRRCIL3','VLIRRCIL4','VLIRRCIL5','VLIRRCIL6','VLIRRCIL7','VLIRRCIL8','VLIRRCIL9','VLIRRCIL10'};
VLIRRCIU={'VLIRRCIU1','VLIRRCIU2','VLIRRCIU3','VLIRRCIU4','VLIRRCIU5','VLIRRCIU6','VLIRRCIU7','VLIRRCIU8','VLIRRCIU9','VLIRRCIU10'};
VLIRRp={'VLIRRp1','VLIRRp2','VLIRRp3','VLIRRp4','VLIRRp5','VLIRRp6','VLIRRp7','VLIRRp8','VLIRRp9','VLIRRp10'};

data.NonVLCase1=data.StudyPop1-data.VLCase1;
n=size(data,1);
data.VLOR1=NaN(n,1);
data.VLORCIL1=NaN(n,1);
data.VLORCIU1=NaN(n,1);
data.VLp1=NaN(n,1);
data.VLRR1=NaN(n,1);
data.VLRRCIL1=NaN(n,1);
data.VLRRCIU1=NaN(n,1);
data.VLRRp1=NaN(n,1);
data.VLIRR1=NaN(n,1);
data.VLIRRCIL1=NaN(n,1);
data.VLIRRCIU1=NaN(n,1);
data.VLIRRp1=NaN(n,1);
for i=2:numel(VLCases)
    % OR
    data.(NonVLCases{i})=data.(StudyPops{i})-data.(VLCases{i});
    data.(VLOR{i})=(data.(VLCases{i}).*data.NonVLCase1)./(data.(NonVLCases{i}).*data.VLCase1);
    stderr=sqrt(1./data.(VLCases{i})+1./data.(NonVLCases{i})+1./data.VLCase1+1./data.NonVLCase1);
    z=log(data.(VLOR{i}))./stderr;
    data.(VLp{i})=2*normcdf(-abs(z));
    data.(VLORCIL{i})=exp(log(data.(VLOR{i}))-zalpha*stderr);
    data.(VLORCIU{i})=exp(log(data.(VLOR{i}))+zalpha*stderr);
    % RR
    data.(VLRR{i})=(data.(VLCases{i}).*data.StudyPop1)./(data.(StudyPops{i}).*data.VLCase1);
    stderrRR=sqrt(1./data.(VLCases{i})-1./data.(StudyPops{i})+1./data.VLCase1-1./data.StudyPop1);
    zRR=log(data.(VLRR{i}))./stderrRR;
    data.(VLRRp{i})=2*normcdf(-abs(zRR));
    data.(VLRRCIL{i})=exp(log(data.(VLRR{i}))-zalpha*stderrRR);
    data.(VLRRCIU{i})=exp(log(data.(VLRR{i}))+zalpha*stderrRR);
    % IRR
    data.(VLIRR{i})=data.(VLIncs{i})./data.VLInc1;
    stderrIRR=sqrt(1./data.(VLCases{i})+1./data.VLCase1);
    zIRR=log(data.(VLIRR{i}))./stderrIRR;
    data.(VLIRRp{i})=2*normcdf(-abs(zIRR));
    data.(VLIRRCIL{i})=exp(log(data.(VLIRR{i}))-zalpha*stderrIRR);
    data.(VLIRRCIU{i})=exp(log(data.(VLIRR{i}))+zalpha*stderrIRR);
end

save('data.mat','data')
