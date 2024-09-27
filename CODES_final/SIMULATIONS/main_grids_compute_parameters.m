clear 
%Compute UR and estimate parameters with grids data from 1880-2022
p1 = genpath('d:\WORKA\Climate\UR_letter\CODES_final\functions');
path(path,p1) %auxiliary functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Working with original data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lectura de datos originales
load DATA_GRID
DATA=DATA_GRID.temp;
years=unique(DATA_GRID.years);
nyears=length(years);
ns=size(DATA,2);

G=NaN(nyears-1,ns);
for i=1:nyears-1
    Y=DATA((i-1)*12+1:i*12,:);
    G(i,:)=nanmean(Y);
end
GRIDS=G(31:end,:); %start in 1880, end in 2022

[T,N]=size(GRIDS);
PARAM.GRIDS.T=T;
PARAM.GRIDS.N=N;

% Select stations present in 1880-2022
GRIDS_1880_2022 = GRIDS(:,all(~isnan(GRIDS)));
G = GRIDS_1880_2022;
[t,n]=size(G);
PARAM.selected_grids.t=t;
PARAM.selected_grids.n=n;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               0.-Test for UR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmax=round(12*(t/100)^0.25);
ADF=NaN(n,1);
UR=NaN(n,1);
for i=1:n
    g=G(:,i);
    [adf,p]=calcula_adf_t_sbic(g,pmax); %ADF with intercep and trend and SBIC
    ADF(i)=adf;
    [h,pValue,stat,cValue,reg]=adftest(g,'model','TS','lags',p); %Compute ADF
    if pValue<=0.1
        UR(i)=0;
    else
        UR(i)=1;
    end
end
disp('The % of grids with UR is:')
sum(UR)/n
PARAM.UR.ADF=ADF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   0.- Looking for SB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model=3: Y{t}=a0+a1*DU+b0*t+b1*DT+e{t} where DU=1(t>TB) and DT=1(t>TB)*(t-TB)
model=3;  
criteria=2; %BIC if criteria=2, AIC if criteria=1
eps=0.15;	%chose a value from {0.01, 0.05, 0.10, 0.15, 0.25}
SB=zeros(n,1);
TB=NaN(n,1);
CV=NaN(n,1);
PY=NaN(n,1);
for i=1:n
    [test,cv,tb]=qfgls_simus(G(:,i),pmax,model, criteria,eps);
    PY(i)=test;
    CV(i)=cv(2);
    if PY(i)>CV(i)
        SB(i)=1;
        TB(i)=tb;
    else
        SB(i)=0;
    end
end
disp('The number of series with SB is:')
sum(SB)/n

PARAM.SB.PY=PY;
PARAM.SB.TB=TB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   1.- Linear trend model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trend=(1:t)';
X=[ones(t,1),trend];
ALPHA=NaN(n,1);
BETA=NaN(n,1);
TBETA=NaN(n,1);
AR=NaN(n,1);
for i=1:n
    g=G(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    ALPHA(i)=b(1);
    BETA(i)=b(2);
    TBETA(i)=t_nw(2);
    [phi, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(res,1,1); %Estimate ar 
    AR(i)=phi(2);
end

PARAM.linear_model_1880_2022.ALPHA=ALPHA;
PARAM.linear_model_1880_2022.BETA=BETA;
PARAM.linear_model_1880_2022.TBETA=TBETA;
PARAM.linear_model_1880_2022.AR=AR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   2.- Cuadratic trend model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[ones(t,1),trend,trend.^2];
ALPHA=NaN(n,1);
BETA1=NaN(n,1);
BETA2=NaN(n,1);
TBETA1=NaN(n,1);
TBETA2=NaN(n,1);
AR=NaN(n,1);
CT2=zeros(n,1);
for i=1:n
    g=G(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    ALPHA(i)=b(1);
    BETA1(i)=b(2);
    BETA2(i)=b(3);
    TBETA1(i)=t_nw(2);
    TBETA2(i)=t_nw(3);
    [phi, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(res,1,1); %Estimate ar 
    AR(i)=phi(2);
    if abs(TBETA2(i))>1.96
        CT2(i)=1;
    end
end

disp('The number of models with signficant cuadratic trend is:')
sum(CT2)/n

PARAM.cuadratic_model_1880_2022.ALPHA=ALPHA;
PARAM.cuadratic_model_1880_2022.BETA1=BETA1;
PARAM.cuadratic_model_1880_2022.BETA2=BETA2;
PARAM.cuadratic_model_1880_2022.TBETA1=TBETA1;
PARAM.cuadratic_model_1880_2022.TBETA2=TBETA2;
PARAM.cuadratic_model_1880_2022.AR=AR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   3.- Cubic trend model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[ones(t,1),trend,trend.^2,trend.^3];
ALPHA=NaN(n,1);
BETA1=NaN(n,1);
BETA2=NaN(n,1);
BETA3=NaN(n,1);
TBETA1=NaN(n,1);
TBETA2=NaN(n,1);
TBETA3=NaN(n,1);
AR=NaN(n,1);
CT3=zeros(n,1);
for i=1:n
    g=G(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    ALPHA(i)=b(1);
    BETA1(i)=b(2);
    BETA2(i)=b(3);
    BETA3(i)=b(4);
    TBETA1(i)=t_nw(2);
    TBETA2(i)=t_nw(3);
    TBETA3(i)=t_nw(4);
    [phi, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(res,1,1); %Estimate ar 
    AR(i)=phi(2);
    if abs(TBETA3(i))>1.96
        CT3(i)=1;
    end
end

disp('The number of models with signficant cubic trend is:')
sum(CT3)/n
PARAM.cubic_model_1880_2022.ALPHA=ALPHA;
PARAM.cubic_model_1880_2022.BETA1=BETA1;
PARAM.cubic_model_1880_2022.BETA2=BETA2;
PARAM.cubic_model_1880_2022.BETA3=BETA3;
PARAM.cubic_model_1880_2022.TBETA1=TBETA1;
PARAM.cubic_model_1880_2022.TBETA2=TBETA2;
PARAM.cubic_model_1880_2022.TBETA3=TBETA3;
PARAM.cubic_model_1880_2022.AR=AR;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   4.- Model with structural break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALPHA1=NaN(n,1);
ALPHA2=NaN(n,1);
BETA1=NaN(n,1);
BETA2=NaN(n,1);
AR=NaN(n,1);
for i=1:n
    g=G(:,i);
    if isnan(TB(i))==0
    tb=TB(i);
    Du=[zeros(tb,1);ones(t-tb,1)];
    Dt=[zeros(tb,1);trend(tb+1:end)-tb];
    X=[ones(t,1),Du,trend,Dt];
    [b,~, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate model with SB
    ALPHA1(i)=b(1);
    ALPHA2(i)=b(2);
    BETA1(i)=b(3);
    BETA2(i)=b(4);
    [phi, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(res,1,1); %Estimate ar 
    AR(i)=phi(2);
    end
end

PARAM.SB_1880_2022.ALPHA1=ALPHA1;
PARAM.SB_1880_2022.ALPHA2=ALPHA2;
PARAM.SB_1880_2022.BETA1=BETA1;
PARAM.SB_1880_2022.BETA2=BETA2;
PARAM.SB_1880_2022.AR=AR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   5.- Linear trend model since 1960
% We work with observations present in the complete sample, 1880-2022, since 1960.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ALPHA=NaN(n,1);
BETA=NaN(n,1);
AR=NaN(n,1);
G2=G(81:t,:);
t2=size(G2,1);
trend=(1:t2)';
X=[ones(t2,1),trend];
for i=1:n
    g=G2(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    ALPHA(i)=b(1);
    BETA(i)=b(2);
    TBETA(i)=t_nw(2);
    [phi, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(res,1,1); %Estimate ar 
    AR(i)=phi(2);
end

PARAM.linear_model_1960_2022.all.ALPHA=ALPHA;
PARAM.linear_model_1960_2022.all.BETA=BETA;
PARAM.linear_model_1960_2022.all.TBETA=TBETA;
PARAM.linear_model_1960_2022.all.AR=AR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      6.- Only observations from 1960
%   We work with observations present in the period 1960-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIDS_1960_2022 = GRIDS(81:end,:);
G3 = GRIDS_1960_2022(:,all(~isnan(GRIDS_1960_2022)));
[t3,n3]=size(G3);
ALPHA=NaN(n3,1);
BETA=NaN(n3,1);
AR=NaN(n3,1);
ADF=NaN(n3,1);
trend=(1:t3)';
X=[ones(t3,1),trend];
for i=1:n3
    g=G3(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    ALPHA(i)=b(1);
    BETA(i)=b(2);
    TBETA(i)=t_nw(2);
    [phi, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(res,1,1); %Estimate ar 
    AR(i)=phi(2);
end

PARAM.linear_model_1960_2022.only.ALPHA=ALPHA;
PARAM.linear_model_1960_2022.only.BETA=BETA;
PARAM.linear_model_1960_2022.only.TBETA=TBETA;
PARAM.linear_model_1960_2022.only.AR=AR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          7.-New observations since 1960
% We work with observations that are new in the period 1960-2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos_1880 = find(all(~isnan(GRIDS))==1);
pos_1960 = find(all(~isnan(GRIDS_1960_2022))==1);
from_1880 = ismember(pos_1960, pos_1880);
mean(BETA(from_1880 == 0));
mean(BETA(from_1880 == 1));

PARAM.linear_model_1960_2022.new.BETA=BETA(from_1880 == 0);
%Notice that is the same that case 5


save('PARAM_GRIDS.mat','PARAM')