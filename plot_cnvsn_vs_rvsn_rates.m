function plot_cnvsn_vs_rvsn_rates(p,d,r,l,pcr,varargin)
%PLOT_CNVSN_VS_RVSN_RATES Scatter plot the estimates of gamma and lambda
pp = table2array(p(:,5:end));
figure;
% DAT studies
h1=loglog(pp(d,1),pp(d,7),'bx',[pp(d,1) pp(d,1)]',[pp(d,end-3) pp(d,end-2)]','b-',[pp(d,2) pp(d,3)]',[pp(d,7) pp(d,7)]','b-'); hold on
% rK39 studies
h2=loglog(pp(r,1),pp(r,7),'kx',[pp(r,1) pp(r,1)]',[pp(r,end-3) pp(r,end-2)]','k-',[pp(r,2) pp(r,3)]',[pp(r,7) pp(r,7)]','k-');
% LST studies
h3=loglog(pp(l,1),pp(l,7),'gx',[pp(l,1) pp(l,1)]',[pp(l,end-3) pp(l,end-2)]','g-',[pp(l,2) pp(l,3)]',[pp(l,7) pp(l,7)]','g-');
% PCR studies
h4=loglog(pp(pcr,1),pp(pcr,7),'cx',[pp(pcr,1) pp(pcr,1)]',[pp(pcr,end-3) pp(pcr,end-2)]','c-',[pp(pcr,2) pp(pcr,3)]',[pp(pcr,7) pp(pcr,7)]','c-');
% Estimates from longitudinal data
% Hasker DAT:
lH=-log(1-252/9873);
[lb,ub]=poiss_CI(252,9873);
lHCI=-log(1-[lb,ub]);
gH=-log(1-216/664);
[lb,ub]=poiss_CI(216,664);
gHCI=-log(1-[lb,ub]);
% Hasker rK39:
lHr=-log(1-145/9873);
[lb,ub]=poiss_CI(145,9873);
lHrCI=-log(1-[lb,ub]);
gHr=-log(1-372/626);
[lb,ub]=poiss_CI(372,626);
gHrCI=-log(1-[lb,ub]);
% Ostyn DAT:
lO=-log(1-375/9034);
[lb,ub]=poiss_CI(375,9034);
lOCI=-log(1-[lb,ub]);
lO2=-log(1-130/8617);
gO=-log(1-318/368);
[lb,ub]=poiss_CI(318,368);
gOCI=-log(1-[lb,ub]);
% Bern rK39:
lB=63.1/1000;
gB=502.1/1000;
h5=loglog(lH,gH,'c.',lHr,gHr,'r.',lO,gO,'m.',lB,gB,'k.',lHCI,[gH gH],'c-',[lH lH],gHCI,'c-',lHrCI,[gHr gHr],'r-',[lHr lHr],gHrCI,'r-',lOCI,[gO gO],'m-',[lO lO],gOCI,'m-','MarkerSize',20);
set(gca,'FontSize',16)
xlabel('\lambda (yr^{-1})')%,'FontSize',24); 
ylabel('\gamma (yr^{-1})')%,'FontSize',24); 
legend([h1(1);h2(1);h3(1);h4(1);h5(1:4)],'DAT','rK39','LST','PCR','Hasker 2013 DAT','Hasker 2013 rK39','Ostyn 2011 DAT','Bern 2007 rK39','Location','NorthWest'); hold off
axis([2e-4 2 4e-3 80])
% axis([0.95*min(pp(:,2)) 1.05*max(pp(:,3)) -Inf Inf])
if nargin==6
    saveas(gcf,varargin{1})
    saveas(gcf,[varargin{1} '.eps'],'epsc')
end