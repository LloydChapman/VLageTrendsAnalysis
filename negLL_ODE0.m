function L = negLL_ODE0(params,data,cens,freq,b0,b1)
% Likelihood function without accounting for averaging over each age group

% Reshape data input from vector to matrix
data=[data(1:4:end),data(2:4:end),data(3:4:end),data(4:4:end)]';

if isempty(b0) && isempty(b1)
    b = params(1:2);
    gamma = params(3);
elseif ~isempty(b0) && isempty(b1)
    b = [b0 params(1)];
    gamma = params(2);
elseif isempty(b0) && ~isempty(b1)
    b = [params(1) b1];
    gamma = params(2);
elseif ~isempty(b0) && ~isempty(b1)
    b = [b0 b1];
    gamma = params(1);
end

a = (data(1,:)+data(2,:)+1)/2; n = data(3,:); k = data(4,:);
sol = ode45(@(t,Y)Rev_Cat_ODE(t,Y,b,gamma),[0 100],0);
p = deval(sol,a);
L_a =  k.*log(p) + (n-k).*log(1-p);
L = -sum(L_a);