function L = negLL_Hasker(pars,data,cens,freq,b0,b1)
% Likelihood function without accounting for averaging over each age group

% Reshape data input from vector to matrix
data=[data(1:6:end),data(2:6:end),data(3:6:end),data(4:6:end),data(5:6:end),data(6:6:end)]';

if isempty(b0) && isempty(b1)
    b = exp(pars(1:2));
    gamma = exp(pars(3));
elseif ~isempty(b0) && isempty(b1)
    b = [b0 exp(pars(1))];
    gamma = exp(pars(2));
elseif isempty(b0) && ~isempty(b1)
    b = [exp(pars(1)) b1];
    gamma = exp(pars(2));
elseif ~isempty(b0) && ~isempty(b1)
    b = [b0 b1];
    gamma = exp(pars(1));
end

a = (data(1,:)+data(2,:)+1)/2; n = data(3,:); k = data(4,:); m = data(5,:); q=data(6,:);
sol = ode45(@(t,Y)Rev_Cat_ODE(t,Y,b,gamma),[0 100],0);
p = deval(sol,a);
L_a =  k.*log(p) + (n-k).*log(1-p) + q.*log(age_dep_FOI(a,b).*(1-p).*m) - age_dep_FOI(a,b).*(1-p).*m;
L = -sum(L_a);
end