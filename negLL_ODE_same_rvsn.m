function [NLL_tot,NLL] = negLL_ODE_same_rvsn(params,data,cens,freq,b0,b1,nagps)
%NEGLL_ODE_SAME_RVSN Negative log-likelihood function for same reversion rate model

m=numel(nagps);
data1=cell(m,1);
% Reshape data input from vector to matrix
for i=1:m
    j=4*sum(nagps(1:i-1));
    k=4*sum(nagps(1:i));
    data1{i}=[data(j+1:4:k),data(j+2:4:k),data(j+3:4:k),data(j+4:4:k)]';
end

if isempty(b0) && isempty(b1)
    b = exp([params(1:m);params(m+1:2*m)]');
elseif ~isempty(b0) && isempty(b1)
    b = [b0*ones(1,m);exp(params(1:m))]';
elseif isempty(b0) && ~isempty(b1)
    b = [exp(params(1:m));b1*ones(1,m)]';
elseif ~isempty(b0) && ~isempty(b1)
    b = [b0*ones(m,1);b1*ones(m,1)];
end
gamma=exp(params(end));

NLL=NaN(m,1);
for i=1:m
    temp = data1{i};
    a = (temp(1,:)+temp(2,:)+1)/2;
    n = temp(3,:);
    k = temp(4,:);
    sol = ode45(@(t,Y)Rev_Cat_ODE(t,Y,b(i,:),gamma),[0 100],0);
    p = deval(sol,a);
    LL_a =  k.*log(p) + (n-k).*log(1-p);
    NLL(i) = -sum(LL_a);
end
NLL_tot=sum(NLL);