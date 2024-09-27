function test=UR3k(y,That,k)

% This file is for the UR test described in Kim and Perron (2009).
% "y" is a column vector of observations.
% "That" is the break date.
% "k" is the number of the lags of the difference of the dependent variable
% to be included in the regression.
% This file computes the t-statistic testing for a unit root
% using Model A3 in Kim and Perron (2009).
% See equation (2) in Kim and Perron (2009) for the description of
% the procedure.

T=size(y,1);
z1=[ones(T,1),(1:T)'];
z2=[zeros(That,1);ones(T-That,1)];
z3=[zeros(That,1);(1:T-That)'];
z=[z1,z2,z3];

yhat=y-z*((z'*z)\(z'*y));

Dy=zeros(T-k-1,k);
for idn=1:k
    Dy(:,idn)=yhat(k+2-idn:T-idn,1)-yhat(k+1-idn:T-idn-1,1);
end
Dt=zeros(T-k-1,k+1);
for ind=0:k
    Dt(That-k+ind,ind+1)=1;
end
yh=yhat(k+2:T,1);
yl=yhat(k+1:T-1,1);
if k==0;
    W=Dt;
else
    W=[Dy,Dt];
end
yh=yh-W*((W'*W)\(W'*yh));
yl=yl-W*((W'*W)\(W'*yl));
   
%Computing t-statistic
a=(yl'*yh)/(yl'*yl);
s2=(yh-a*yl)'*(yh-a*yl)/(T-3-3*k);


test=(a-1)/sqrt(s2/(yl'*yl));
