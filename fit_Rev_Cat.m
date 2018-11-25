function [par,parCI,NLL] = fit_Rev_Cat(data,str,doPlots,b0,b1,test)
%FIT_REV_CAT Fit reversible catalytic model

% If no arguments are specified use Hasker 2013 data
if nargin == 0
    data = [2 10 20 30 40 50 60 70 ; 
            9 19 29 39 49 59 69 90 ;
            3858 2802 1565 1459 1021 812 767 321 ;
            100 126 92 123 101 87 97 51];
    str = 'Hasker et al, 2013';
    doPlots = true;
    b0 = [];
    b1 = [];
end

par0 = [0.01 0.001 0.1];
if isempty(b0) && isempty(b1)
    logpar0 = log(par0);
elseif ~isempty(b0) && isempty(b1)
    logpar0 = log(par0(2:3));
elseif isempty(b0) && ~isempty(b1)
    logpar0 = log(par0([1,3]));
elseif ~isempty(b0) && ~isempty(b1)
    logpar0 = log(par0(3));
end

if ~(strcmp(str,'Nandy et al, 1987') && ~isempty(b1))
    [logpar, logparCI] = mle(data(:),'nloglf',@(params,data,cens,freq)negLL_ODE(params,data,cens,freq,b0,b1),'start',logpar0);
    par = exp(logpar);
    parCI = exp(logparCI);
    NLL = negLL_ODE(logpar,data(:),[],[],b0,b1);
else % if fitting const. rate model to Nandy data don't do log transformation of parameters as 95% CI overlaps gamma=0 
    [par,parCI] = mle(data(:),'nloglf',@(params,data,cens,freq)negLL_ODE0(params,data,cens,freq,[],0),'start',par0([1,3]),'optimfun','fmincon'); 
    NLL = negLL_ODE0(par,data(:),[],[],[],0);
    parCI(parCI<0)=0;
end

if doPlots
    figure;
    h=plot_model_fit(data,par,b0,b1);
    legend(h,'Data','Model')
    title(str)
    saveas(gcf,[char(str) test '.eps'],'epsc')
    
    if isempty(b0) && isempty(b1)
        plot_LL_ODE(data,str,par,parCI,NLL)
    elseif isempty(b0) && b1==0
        plot_LL(data,str,par,parCI,NLL)
    end
end

end