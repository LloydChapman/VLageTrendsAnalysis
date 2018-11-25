function plot_LL(data,str,par,parCI,NLL)
%PLOT_LL Plot log-likelihood surface for different values of lambda (b0) 
%and gamma for age-independent infection rate model
%   The 2-d 95% profile confidence interval and the 95% confidence interval
%   based on a Normal approximation to log-likelihood surface at the MLE
%   using the Hessian are plotted.
    lvec = linspace(parCI(1,1)/2,2*parCI(2,1));
    gvec = linspace(parCI(1,2)/2,2*parCI(2,2));
%     lvec = linspace(parCI(1,1)/10,10*parCI(2,1));
%     gvec = linspace(parCI(1,2)/10,10*parCI(2,2));
    LLmat = zeros(numel(gvec),numel(lvec));
    for i = 1:numel(lvec)
        for j = 1:numel(gvec)
            LLmat(j,i) = -negLL_ODE(log([lvec(i) gvec(j)]),data(:),[],[],[],0);
        end
    end
    figure;
    surf(lvec,gvec,LLmat,'EdgeColor','none'); hold on
    plot3(par(1),par(2),-NLL,'r.',[parCI(1,1) parCI(1,1)],gvec([1,end]),-[NLL NLL],'r-',[parCI(2,1) parCI(2,1)],gvec([1,end]),-[NLL NLL],'r-',lvec([1,end]),[parCI(1,2) parCI(1,2)],-[NLL NLL],'r-',lvec([1,end]),[parCI(2,2) parCI(2,2)],-[NLL NLL],'r-','LineWidth',2,'MarkerSize',20); hold on
    v = -NLL-chi2inv(0.95,1)/2;
    contour3(lvec,gvec,LLmat,[v v],'w','LineWidth',2)
    hold off
%     xlim(parCI(:,1))
%     ylim(parCI(:,2))
    xlabel('\lambda')
    ylabel('\gamma')
    zlabel('LL')
    title(str)
    view([0 90])
end