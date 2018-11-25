function p=solve_Rev_Cat_ODE(a,b,gamma)
sol = ode45(@(t,Y)Rev_Cat_ODE(t,Y,b,gamma),[0 100],0); 
p = deval(sol,a);