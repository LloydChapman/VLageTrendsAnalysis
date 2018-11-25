function [par,parCI,NLL] = fit_Rev_Cat_Hasker(data,str,doPlots,b0,b1)
%FIT_REV_CAT_HASKER Fit reversible catalytic model to Hasker 2013
%seropositivity and seroconversion data

% If no arguments are specified, use Hasker 2013 data 
if nargin == 0
    data = [2 10 20 30 40 50 60 70 ; 
            9 19 29 39 49 59 69 90 ;
            3858 2802 1565 1459 1021 812 767 321 ;
            100 126 92 123 101 87 97 51;
            3415 2146 1093 1084 765 602 555 213 ;
            65 59 54 44 45 38 31 15];
    str = 'Hasker et al, 2013';
    doPlots = true;
    b0 = [];
    b1 = [];
end

par0 = [0.01 0.001 0.5];
if isempty(b0) && isempty(b1)
    logpar0 = log(par0);
elseif ~isempty(b0) && isempty(b1)
    logpar0 = log(par0(2:3));
elseif isempty(b0) && ~isempty(b1)
    logpar0 = log(par0([1,3]));
elseif ~isempty(b0) && ~isempty(b1)
    logpar0 = log(par0(3));
end
[logpar, logparCI] = mle(data(:),'nloglf',@(pars,data,cens,freq)negLL_Hasker(pars,data,cens,freq,b0,b1),'start',logpar0);
par = exp(logpar);
parCI = exp(logparCI);
NLL = negLL_Hasker(logpar,data(:),[],[],b0,b1);

a_m = (data(1,:)+data(2,:)+1)/2;
p_data = data(4,:)./data(3,:);
% l_data = data(6,:)./data(5,:);
l_data = -log(1-data(6,:)./data(5,:)); %N.B. check the point estimates should be logged to convert them to continuous rates
if isempty(b0) && isempty(b1)
    b = par(1:2);
    gamma = par(3);
elseif ~isempty(b0) && isempty(b1)
    b = [b0 par(1)];
    gamma = par(2);
elseif isempty(b0) && ~isempty(b1)
    b = [par(1) b1];
    gamma = par(2);
elseif ~isempty(b0) && ~isempty(b1)
    b = [b0 b1];
    gamma = par(1);
end
sol = ode45(@(t,Y)Rev_Cat_ODE(t,Y,b,gamma),[0 100],0);
a = 0:0.1:100;
p = deval(sol,a);
l = age_dep_FOI(a,b);

if doPlots
    figure;
    plot(a_m, p_data,'kx', a, p, 'LineWidth',1.5,'MarkerSize',8); hold on
    %plot the CI for each data point...
    for ia = 1 : length(a_m)
        pd = fitdist(data(4,ia),'binomial','NTrials',data(3,ia));
        ci = paramci(pd,'Parameter','p');
        plot([a_m(ia) a_m(ia)],[ci(1) ci(2)],'k-','LineWidth',1.5);
    end
    set(gca,'FontSize',16)
    xlabel('Age (yrs)');
    ylabel('Proportion DAT positive')
    legend('Data', 'Model');
    title(str)
    saveas(gcf,['HaskerDATprev' num2str(b1) '.eps'],'epsc')
    figure;
    plot(a_m, l_data,'kx', a, l, 'LineWidth',1.5,'MarkerSize',8); hold on
    %plot the CI for each data point...
    for ia = 1 : length(a_m)
        [sCILB,sCIUB] = poiss_CI(data(6,ia),data(5,ia),0.05);
        plot([a_m(ia) a_m(ia)],-log(1-[sCILB sCIUB]),'k-','LineWidth',1.5);
    end
    set(gca,'FontSize',16)
    xlabel('Age (yrs)');
    ylabel('Seroconversion rate (yr^{-1})')
    legend('Data', 'Model');
    title(str)
    saveas(gcf,['HaskerSerocnvsn' num2str(b1) '.eps'],'epsc')
end

end