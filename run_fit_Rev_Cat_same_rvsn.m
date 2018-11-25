function [p,p1,p2]=run_fit_Rev_Cat_same_rvsn(data,DATidx,rK39idx,LSTidx,PCRidx,b0,b1,fname)
%RUN_FIT_REV_CAT_SAME_RVSN Run fitting of same reversion rate reversible catalytic model

% If no arguments supplied, load latest saved version of database and use 
% studies included in review
if nargin==0
    load data
    DATidx=find((strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (DAT)'))|strcmp(data.Author,'Koirala et al, 2004')|(strcmp(data.Author,'Ostyn et al, 2015')&strcmp(data.Type,'PREVALENCE'))|strcmp(data.Author,'Rijal et al, 2010')|(strcmp(data.Author,'Singh et al, 2010')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (DAT)'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (DAT)')));
    rK39idx=find(strcmp(data.Author,'Bern et al, 2007')|(strcmp(data.Author,'Hasker et al, 2013')&strcmp(data.Type,'PREVALENCE (rK39)'))|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (rK39)')));
    LSTidx=find(strcmp(data.Author,'Bern et al, 2006')|(strcmp(data.Author,'Nandy et al, 1987')&strcmp(data.Type,'PREVALENCE'))|(strcmp(data.Author,'Patil et al, 2013')&strcmp(data.Type,'PREVALENCE (LST)'))|(strcmp(data.Author,'Schenkel et al, 2006')&strcmp(data.Type,'PREVALENCE (LST)')));
    PCRidx=find(strcmp(data.Author,'Kaushal et al, 2017')|(strcmp(data.Author,'Topno et al, 2010')&strcmp(data.Type,'PREVALENCE (PCR)')));
    b0=[];
    b1=0;
    fname='RsltsAgeIndepFOISameRvsn';
end

%% Fit to all data with same reversion rate, study-specific conversion rates
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
y=[];
nagps=NaN(1,n);
str=cell(1,n);
count=zeros(ntests,1);
for j=1:n
    i=idx(j);
    dataset=['data' num2str(i) '.mat'];
    if exist(dataset,'file')
        load(dataset)
        str{j}=data.Author(i);
        y=[y x];
        nagps(j)=size(x,2);
        if ismember(i,DATidx)
            count(1)=count(1)+1;
            testdata{1} = [testdata{1} x];
        elseif ismember(i,rK39idx)
            count(2)=count(2)+1;
            testdata{2} = [testdata{2} x];
        elseif ismember(i,LSTidx)
            count(3)=count(3)+1;
            testdata{3} = [testdata{3} x];
        elseif ismember(i,PCRidx)
            count(4)=count(4)+1;
            testdata{4} = [testdata{4} x];
        end
    end
end
[MLE,CI,NLL]=fit_Rev_Cat_same_rvsn(y,str,false,b0,b1,nagps);
p=put_in_table(p,1:n,MLE,CI,NLL,b0,b1,n);

if isempty(b0) && ~isempty(b1) && b1==0
    p=add_total(p,n+1);
elseif isempty(b0) && isempty(b1)
    p=add_total(p,2*n+1);
end

k=0;
nagps1=cell(ntests,1);
nagps2=NaN(ntests,1);
for j=1:ntests
    nagps1{j}=nagps(k+1:k+count(j));
    k=k+count(j);
    nagps2(j)=sum(nagps1{j});
end

if isempty(b0) && ~isempty(b1) && b1==0
% plot_cnvsn_vs_rvsn_rates(p,d,r,l,pcr)
end

%% Fit to all data with test-specific reversion rates, study-specific conversion rate
p1=p(1:end-1,:);
k=0;
for j=1:ntests
    nstds=numel(nagps1{j});
    [MLE,CI,NLL]=fit_Rev_Cat_same_rvsn(testdata{j},str(k+1:k+nstds),false,b0,b1,nagps1{j});
    p1=put_in_table(p1,k+1:k+nstds,MLE,CI,NLL,b0,b1,nstds);
    k=k+nstds;
end

if isempty(b0) && ~isempty(b1) && b1==0
    p1=add_total(p1,n+ntests);
elseif isempty(b0) && isempty(b1)
    p1=add_total(p1,2*n+ntests);
end

%% Fit to all data with same reversion rate, test-specific conversion rate
Author=cellfun(@(x)[x ' studies'],tests,'UniformOutput',false);
p2=table(Author);
p2.Location=cell(ntests,1);
p2.Year=NaN(ntests,1);
p2.Test=tests;
p2.b0=NaN(ntests,1);
p2.b0LB=NaN(ntests,1);
p2.b0UB=NaN(ntests,1);
p2.b1=NaN(ntests,1);
p2.b1LB=NaN(ntests,1);
p2.b1UB=NaN(ntests,1);
p2.gamma=NaN(ntests,1);
p2.gammaLB=NaN(ntests,1);
p2.gammaUB=NaN(ntests,1);
p2.NLL=NaN(ntests,1);
p2.AIC=NaN(ntests,1);
[MLE,CI,NLL]=fit_Rev_Cat_same_rvsn(y,tests,false,b0,b1,nagps2);
p2=put_in_table(p2,1:ntests,MLE,CI,NLL,b0,b1,ntests);

if isempty(b0) && ~isempty(b1) && b1==0
    p2=add_total(p2,ntests+1);
elseif isempty(b0) && isempty(b1)
    p2=add_total(p2,2*ntests+1);
end

writetable(vertcat(p2,p,p1),[fname '.csv'],'QuoteStrings',true)
print_par_ests(vertcat(p2,p,p1),fname)

end