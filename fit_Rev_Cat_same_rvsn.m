function [par,parCI,NLL] = fit_Rev_Cat_same_rvsn(data,str,doPlots,b0,b1,nagps)
%FIT_REV_CAT_SAME_RVSN Fit same reversion rate reversible catalytic model

% If no arguments are specified, use Hasker 2013 data
if nargin == 0
    data = [2 10 20 30 40 50 60 70 ; 
            9 19 29 39 49 59 69 90 ;
            3858 2802 1565 1459 1021 812 767 321 ;
            100 126 92 123 101 87 97 51];
    str = {'Hasker et al, 2013'};
    doPlots = true;
    b0 = [];
    b1 = 0;
    nagps = 8;
end

m = numel(nagps);
par0 = [0.01*ones(1,m) 0.001*ones(1,m) 0.5];
if isempty(b0) && isempty(b1)
    logpar0 = log(par0);
elseif ~isempty(b0) && isempty(b1)
    logpar0 = log(par0(m+1:end));
elseif isempty(b0) && ~isempty(b1)
    logpar0 = log(par0([1:m,end]));
elseif ~isempty(b0) && ~isempty(b1)
    logpar0 = log(par0(end));
end

[logpar, logparCI] = mle(data(:),'nloglf',@(params,data,cens,freq)negLL_ODE_same_rvsn(params,data,cens,freq,b0,b1,nagps),'start',logpar0);
par = exp(logpar);
parCI = exp(logparCI);
[NLL_tot,NLL] = negLL_ODE_same_rvsn(logpar,data(:),[],[],b0,b1,nagps);

if doPlots
    k=0;
    for i=1:m
        datai=data(:,k+1:k+nagps(i));
        if isempty(b0) && isempty(b1)
            pari = par([i,m+i,end]);
        elseif ~isempty(b0) && isempty(b1)
            pari = par([i,end]);
        elseif isempty(b0) && ~isempty(b1)
            pari = par([i,end]);
        elseif ~isempty(b0) && ~isempty(b1)
            pari = par(1);
        end
        figure;
        h=plot_model_fit(datai,pari,b0,b1);
        legend(h,'Data','Model')
        title(str{i})
        k=k+nagps(i);
    end
end

end