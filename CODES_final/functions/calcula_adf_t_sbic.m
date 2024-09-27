%May 2008
%THIS PROGRAM ESTIMATES DFA
%LAG SELECTED WITH SBIC

function[dfa,p]=calcula_adf_t_sbic(y,pmax)
%pmax=round(12*(t/100)^0.25);
for i=0:1:pmax
[beta,sigma,res]=ols_adf_t(y,i);    
T=length(res);
sigma2_res=res'*res/T;
sbic(i+1)=log(sigma2_res)+(i+3)*log(T)/T;
end
[sbic p]=min(sbic); %p selected
p=p-1; %correction for the case lag=0
[beta,sigma,res]=ols_adf_t(y,p);
dfa=beta(3)/sigma(3);

%This function estimates ADF regresion with trend
function [beta,sigma,res]=ols_adf_t(y,p)
n=length(y);
    for j=1:n-1
    dy(j)=y(j+1)-y(j);
    end
x=zeros(n-1-p,p);
%let's now fill the columns pn x...
for j=1:p
    x(:,j)=dy(p-j+1:n-1-j);
end
trend=1:1:n-1-p;
x=[ones(n-1-p,1) trend' y(p+1:n-1) x];
[t k]=size(x);
beta=ols(dy(p+1:n-1)',x);
res=dy(p+1:n-1)'-x*beta;
sigma2=inv(x'*x)*(res'*res)/(t-k);
sigma=diag(sigma2).^0.5;

%this function computes OLS estimates using standard techniques

function[beta,beta1,ste,r2,resi]=ols(y,X)

%beta1=(X'*X)\(X'*y);
%a more sophisticated calculation:
[Q, R]=qr(X,0);
beta= R\(Q'*y);
%we continue the program using beta:
%resi=y-X*beta;
%standard errors:
%[n k]=size(X);
%varres=resi'*resi/(n-k);
%varbeta=varres*inv(X'*X);
%ste=diag(varbeta).^0.5;
%r2=1-resi'*resi/((y-mean(y))'*(y-mean(y)));
%Second part:
%compute hi.
%P=X*inv(X'*X)*X';
%h=diag(P);
%v=h./(1-h);
%sum(v);
%plot(v)








