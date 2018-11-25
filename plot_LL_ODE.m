function plot_LL_ODE(data,str,par,parCI,NLL)
%PLOT_LL_ODE Plot log-likelihood surface for each pair of parameters in the
%age-dependent infection rate model
%   The 2-d 95% profile confidence interval and the 95% confidence interval
%   based on a Normal approximation to log-likelihood surface at the MLE
%   using the Hessian are plotted.
parname={'b_0','b_1','\gamma'};
for k=1:3
    m=find(~ismember(1:3,k));
    par1vec = linspace(par(m(1))/2,2*par(m(1)),20);
    par2vec = linspace(par(m(2))/2,2*par(m(2)),20);
    LLmat = zeros(numel(par1vec),numel(par2vec));
    tmp = NaN(1,3);
    for i = 1:numel(par1vec)
        for j = 1:numel(par2vec)
            tmp(m(1)) = par1vec(i);
            tmp(m(2)) = par2vec(j);
            tmp(k) = par(k);
%             LLmat(j,i) = -negLL_ODE([log(tmp(1)) tmp(2) log(tmp(3))],data(:));
            LLmat(j,i) = -negLL_ODE(log(tmp),data(:),[],[],[],[]);
        end
    end
    figure;
    surf(par1vec,par2vec,LLmat,'EdgeColor','none'); hold on
    plot3(par(m(1)),par(m(2)),-NLL,'r.','MarkerSize',30); 
    v = -NLL-chi2inv(0.95,1);
    contour3(par1vec,par2vec,LLmat,[v v],'w','LineWidth',2)
    hold off
%     xlim(parCI(:,1))
%     ylim(parCI(:,2))
    xlabel(parname{m(1)})
    ylabel(parname{m(2)})
    zlabel('LL')
    title(str)
end
end