function dpdt=Rev_Cat_ODE(t,p,b,gamma)
dpdt=age_dep_FOI(t,b)*(1-p)-gamma*p;