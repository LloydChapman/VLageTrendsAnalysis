function p=put_in_table(p,j,MLE,CI,NLL,b0,b1,m)
if isempty(b0) && isempty(b1)
    p.b0(j)=MLE(1:m);
    p.b0LB(j)=CI(1,1:m);
    p.b0UB(j)=CI(2,1:m);
    p.b1(j)=MLE(m+1:2*m);
    p.b1LB(j)=CI(1,m+1:2*m);
    p.b1UB(j)=CI(2,m+1:2*m);
    p.gamma(j)=MLE(end);
    p.gammaLB(j)=CI(1,end);
    p.gammaUB(j)=CI(2,end);
    if m==1
        p.AIC(j)=AIC(NLL,2*m+1);
    else
        p.AIC(j)=NaN(m,1);
    end
elseif ~isempty(b0) && isempty(b1)
    p.b0(j)=b0*ones(m,1);
    p.b1(j)=MLE(1:m);
    p.b1LB(j)=CI(1,1:m);
    p.b1UB(j)=CI(2,1:m);
    p.gamma(j)=MLE(m+1);
    p.gammaLB(j)=CI(1,m+1);
    p.gammaUB(j)=CI(2,m+1);
    if m==1
        p.AIC(j)=AIC(NLL,m+1);
    else
        p.AIC(j)=NaN(m,1);
    end
elseif isempty(b0) && ~isempty(b1)
    p.b0(j)=MLE(1:m);
    p.b0LB(j)=CI(1,1:m);
    p.b0UB(j)=CI(2,1:m);
    p.b1(j)=b1*ones(m,1);
    p.gamma(j)=MLE(m+1);
    p.gammaLB(j)=CI(1,m+1);
    p.gammaUB(j)=CI(2,m+1);
    if m==1
        p.AIC(j)=AIC(NLL,m+1);
    else
        p.AIC(j)=NaN(m,1);
    end
elseif ~isempty(b0) && ~isempty(b1)
    p.b0(j)=b0*ones(m,1);
    p.b1(j)=b1*ones(m,1);
    p.gamma(j)=MLE(1);
    p.gammaLB(j)=CI(1,1);
    p.gammaUB(j)=CI(2,1);
    if m==1
        p.AIC(j)=AIC(NLL,1);
    else
        p.AIC(j)=NaN(m,1);
    end
end
p.NLL(j)=NLL;
