function khat=MIC(y,m,option,Tb)
% Modified Information Criterion for the AO models.

% m 0: DF with an intercept and slope 1: crash, 2: changing growth, 3: mixed model

% option 1: MAIC, 2:MBIC, 3:AIC, 4:BIC

% Tb is the break date. When m=0, Tb can be any positive integer between 0
% and the sample size.

% The output "khat" stands for the number of lags of first difference, NOT AR order.

T=size(y,1);
kmax=10;%round(10*(T/100)^(1/4));

switch option
    case 1
        CT=2;
    case 2
        CT=log(T);
    case 3
        CT=2;
    case 4
        CT=log(T-kmax);
end

z1=[ones(T,1),(1:T)'];
z2=[zeros(Tb,1);ones(T-Tb,1)];
z3=[zeros(Tb,1);(1:T-Tb)'];

switch m
    case 0
        z=z1;
    case 1
        z=[z1,z2];
    case 2
        z=[z1,z3];
    case 3
        z=[z1,z2,z3];
end

%Detrending the given series
yhat=y-z*((z'*z)\(z'*y));

%Regressand
dyhat=yhat(kmax+2:T,1)-yhat(kmax+1:T-1,1);

%Regressors
yhatl=yhat(kmax+1:T-1,1);
Rg=zeros(T-kmax-1,kmax);
for idn=1:kmax
    Rg(:,idn)=yhat(kmax+2-idn:T-idn,1)-yhat(kmax+1-idn:T-1-idn,1);
end

S=zeros(kmax+1,1);
for j=0:kmax
    X=[yhatl,Rg(:,1:j)];
    bhat=(X'*X)\(X'*dyhat);
    b0=bhat(1,1);
    sig2=(dyhat-X*bhat)'*(dyhat-X*bhat)/size(dyhat,1);
    switch option
        case 1
            tau=(b0^2)*yhatl'*yhatl/(T^2)/sig2;
            S(j+1,1)=log(sig2)+CT*(tau+j)/T;
        case 2
            tau=(b0^2)*yhatl'*yhatl/(T^2)/sig2;
            S(j+1,1)=log(sig2)+CT*(tau+j)/T;
        case 3
            S(j+1,1)=log(sig2)+CT*j/(T-kmax);
        case 4
            S(j+1,1)=log(sig2)+CT*j/(T-kmax);
    end
end

[MinS,khat]=min(S);
khat=khat-1;