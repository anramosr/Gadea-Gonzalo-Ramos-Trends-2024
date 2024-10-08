function [That]=breakestimate2(y,m,eps)

% This file is for the estimation of the break date.
% For the precise description of the procedure,
% see Kim and Perron (2009) equation (5).
% "y" is a column vector of observations.
% "m" is for an option. It should be 1, 2 or 3.
% m=1 is for Model A1, m=2 for Model A2, and m=3 for Model A3
% in Kim and Perron (2009).
% eps is the % of observations at the beginning and end of the sample that are discarded to look for breaks.

% The output "That" is the estimated break date.
[T,k]=size(y);
k=round(eps*T); %remove observations at the begining and the end of the sample

z1=[ones(T,1),(1:T)'];

SSR=zeros(T,1);
for Tb=k+1:T-k
    switch m
        case 1
            z2search=[zeros(Tb,1);ones(T-Tb,1)];
            zsearch=[z1,z2search];
        case 2
            z3search=[zeros(Tb,1);(1:T-Tb)'];
            zsearch=[z1,z3search];
        case 3
            z2search=[zeros(Tb,1);ones(T-Tb,1)];
            z3search=[zeros(Tb,1);(1:T-Tb)'];
            zsearch=[z1,z2search,z3search];
    end
    
    ysearch=y-zsearch*((zsearch'*zsearch)\(zsearch'*y));
    ssr=ysearch'*ysearch;
    SSR(Tb,1)=ssr;
end
SSR=SSR(k+1:T-k,1);
[MinSSR,That]=min(SSR);
That=That+k; % The estimate of break date