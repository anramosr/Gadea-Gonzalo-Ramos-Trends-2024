function[beta, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(y,p,constant)
t=length(y);
x=zeros(t-p,p);
%let's now fill the columns in x...
for j=1:p
    x(:,j)=y(p-j+1:t-j);
end

if constant==1
x=[ones(t-p,1) x];
end

[beta,t_nw, se_nw, res, r2, varres, varBhat, y_hat]=ols_hac(y(p+1:t),x);