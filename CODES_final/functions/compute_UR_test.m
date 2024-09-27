function [h_adf,pValue_adf,h_kp]=compute_UR_test(y)


%%%%%%%%%%%%%%ADF_UR_test_trend%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=length(y);
pmax=round(12*(T/100)^0.25);
[adf,p]=calcula_adf_t_sbic(y,pmax); %ADF with intercep and trend and SBIC
[h_adf,pValue_adf,stat,cValue,reg]=adftest(y,'model','TS','lags',p, 'Alpha', 0.01);

%%%%%%%%%%%%%%%%UR with breaks, Kim-Perron%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps=0.15; %eps is the % of observations at the beginning and end of the sample
% %that are discarded to look for breaks. In the original code only two
% observacions are discarted

Dt=breakestimate2(y,3,eps); %estimate break, model 3
khat=MIC(y,0,2,Dt); %select lag with MBIC and model with intercept and slope
%if you use the original you will probably have to make some adjustments between
% the number of lags and the position of the break.
%if khat>=T-Dt khat=T-Dt-1; end
%if khat>=Dt khat=Dt-1; end

test=UR3k(y,Dt,khat); %Unit root test
lambda=round(Dt/T,1); sig=1;
cv=TableIVB_Perron_1989 (lambda,sig);
if abs(test)>abs(cv) h_kp=1; else h_kp=0; end
