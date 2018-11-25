function [CIL,CIU]=PoissCI(N,PT,alpha)
if nargin==2
    alpha=0.05;
end
CIL=chi2inv(alpha/2,2*N)./(2*PT);
CIU=chi2inv(1-alpha/2,2*(N+1))./(2*PT);
end