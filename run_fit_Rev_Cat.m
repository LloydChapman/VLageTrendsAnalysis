function [p,p1,p2]=run_fit_Rev_Cat(data,DATidx,rK39idx,LSTidx,PCRidx,b0,b1,fname)
%RUN_FIT_REV_CAT Run fitting of reversible catalytic model

% If no arguments supplied, load latest saved version of database and only 
% use studies included in review for which the reversible catalytic model
% fits
if nargin==0
    load data
    DATidx=find((strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (DAT)'))|strcmp(data.Author,'Koirala et al, 2004')|(strcmp(data.Author,'Ostyn et al, 2015')&strcmp(data.Type,'PREVALENCE'))|strcmp(data.Author,'Rijal et al, 2010')|(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (DAT)'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (DAT)')));
    rK39idx=find(strcmp(data.Author,'Bern et al, 2007')|(strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (rK39)'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (rK39)')));
    LSTidx=find(strcmp(data.Author,'Bern et al, 2006')|(strcmp(data.Author,'Nandy et al, 1987')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Patil et al, 2013')&strcmp(data.Type,'PREVALENCE (LST)'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (LST)')));
    PCRidx=find(strcmp(data.Author,'Kaushal et al, 2017')|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (PCR)')));
    b0=[];
    b1=0;
    fname='RsltsAgeIndepFOI';
end

%% Fit separately to data from each study
idx=[DATidx;rK39idx;LSTidx;PCRidx];
n=numel(idx);
p=table(data.Author(idx));
p.Properties.VariableNames{'Var1'}='Author';
p.Location=cellfun(@(x,y)[x ', ' y],data.State(idx),data.Country(idx),'UniformOutput',false);
p.Year=data.StartYear(idx);
p.Test=cell(n,1);
d=ismember(idx,DATidx);
r=ismember(idx,rK39idx);
l=ismember(idx,LSTidx);
pcr=ismember(idx,PCRidx);
p.Test(d)={'DAT'};
p.Test(r)={'rK39'};
p.Test(l)={'LST'};
p.Test(pcr)={'PCR'};
p.b0=NaN(n,1);
p.b0LB=NaN(n,1);
p.b0UB=NaN(n,1);
p.b1=NaN(n,1);
p.b1LB=NaN(n,1);
p.b1UB=NaN(n,1);
p.gamma=NaN(n,1);
p.gammaLB=NaN(n,1);
p.gammaUB=NaN(n,1);
p.NLL=NaN(n,1);
p.AIC=NaN(n,1);
tests={'DAT';'rK39';'LST';'PCR'};
ntests=numel(tests);
testdata=cell(ntests,1);
alldata=[];
c=1;
for j=1:n
    i=idx(j);
    dataset=['data' num2str(i) '.mat'];
    if exist(dataset,'file')
        load(dataset)
        %add the data to the list for an overall fit for each test
        if ismember(i,DATidx)
            testdata{1} = [testdata{1} x];
            test = 'DAT';
        elseif ismember(i,rK39idx)
            testdata{2} = [testdata{2} x];
            test = 'rK39';
        elseif ismember(i,LSTidx)
            testdata{3} = [testdata{3} x];
            test = 'LST';
        elseif ismember(i,PCRidx)
            testdata{4} = [testdata{4} x];
            test = 'PCR';
        end
        %add the data to the list for an overall fit
        alldata = [alldata x];
        str=data.Author(i);
        if exist([char(str) test '.eps'],'file')
            test=[test '_' num2str(c)];
            c=c+1;
        end
        [MLE,CI,NLL]=fit_Rev_Cat(x,str,false,b0,b1,test);
        p=put_in_table(p,j,MLE,CI,NLL,b0,b1,1);
    end
end

if isempty(b0) && ~isempty(b1) && b1==0
    p=add_total(p,2*n);
    plot_cnvsn_vs_rvsn_rates(p,d,r,l,pcr,'lambda_gamma_plot_age_indep')
elseif isempty(b0) && isempty(b1)
    p=add_total(p,3*n);
    p1=p;
    a1=20;
    ib1=(p1.b1>1e-12);
    p1.b0(ib1)=p1.b0(ib1)+p1.b1(ib1)*a1;
    p1.b0LB(ib1)=p1.b0LB(ib1)+p1.b1LB(ib1)*a1;
    p1.b0UB(ib1)=p1.b1UB(ib1)+p1.b1UB(ib1)*a1;
    plot_cnvsn_vs_rvsn_rates(p1,d,r,l,pcr,'lambda_gamma_plot_age_dep');
end

%% Fit to data aggregated for each test
Author=cellfun(@(x)[x ' studies'],tests,'UniformOutput',false);
p1=table(Author);
p1.Location=cell(ntests,1);
p1.Year=NaN(ntests,1);
p1.Test=tests;
p1.b0=NaN(ntests,1);
p1.b0LB=NaN(ntests,1);
p1.b0UB=NaN(ntests,1);
p1.b1=NaN(ntests,1);
p1.b1LB=NaN(ntests,1);
p1.b1UB=NaN(ntests,1);
p1.gamma=NaN(ntests,1);
p1.gammaLB=NaN(ntests,1);
p1.gammaUB=NaN(ntests,1);
p1.NLL=NaN(ntests,1);
p1.AIC=NaN(ntests,1);
for j=1:ntests
    [MLE,CI,NLL]=fit_Rev_Cat(testdata{j},['All ' tests{j} ' data'],false,b0,b1,'');
    p1=put_in_table(p1,j,MLE,CI,NLL,b0,b1,1);
end

if isempty(b0) && ~isempty(b1) && b1==0
    p1=add_total(p1,2*ntests);
elseif isempty(b0) && isempty(b1)
    p1=add_total(p1,3*ntests);
end

%% Fit to all data
p2=table({'All data'});
p2.Properties.VariableNames{'Var1'}='Author';
p2.Location={''};
p2.Year=NaN;
p2.Test='All';
p2.b0=NaN;
p2.b0LB=NaN;
p2.b0UB=NaN;
p2.b1=NaN;
p2.b1LB=NaN;
p2.b1UB=NaN;
p2.gamma=NaN;
p2.gammaLB=NaN;
p2.gammaUB=NaN;
p2.NLL=NaN;
p2.AIC=NaN;
[MLE,CI,NLL]=fit_Rev_Cat(alldata,'All data',false,b0,b1,'');
p2=put_in_table(p2,1,MLE,CI,NLL,b0,b1,1);

writetable(vertcat(p2,p1,p),[fname '.csv'],'QuoteStrings',true)
print_par_ests(vertcat(p2,p1,p),fname)

end
