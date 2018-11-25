function p=add_total(p,np)
p.Author{end+1}='Total';
p.Location{end}='';
p.Year(end)=NaN;
p.Test{end}='';
p.b0(end)=NaN;
p.b0LB(end)=NaN;
p.b0UB(end)=NaN;
p.b1(end)=NaN;
p.b1LB(end)=NaN;
p.b1UB(end)=NaN;
p.gamma(end)=NaN;
p.gammaLB(end)=NaN;
p.gammaUB(end)=NaN;
p.NLL(end)=sum(p.NLL);
p.AIC(end)=AIC(p.NLL(end),np);
    
    