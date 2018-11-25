function h = plot_model_fit(data,par,b0,b1)
a_m = (data(1,:)+data(2,:)+1)/2;
p_data = data(4,:)./data(3,:);
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
a = 0:100; 
p = solve_Rev_Cat_ODE(a,b,gamma);

h=plot(a_m, p_data,'kx', a, p,'LineWidth',1.5,'MarkerSize',8); hold on
%plot the CI for each data point...
for ia = 1 : length(a_m)
    pd = fitdist(data(4,ia),'binomial','NTrials',data(3,ia));
    ci = paramci(pd,'Parameter','p');
    plot([a_m(ia) a_m(ia)],[ci(1) ci(2)],'k-');
end
yl = ylim;
ylim([0 yl(2)]);
set(gca,'FontSize',16)
xlabel('Age (yrs)','FontSize',16)
ylabel('Proportion positive','FontSize',16)