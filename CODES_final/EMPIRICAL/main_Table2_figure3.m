clear 
%Compute UR and estimate parameters with grids data from 1880-2022
p1 = genpath('c:\WORKA\Climate\UR_letter\CODES_final\functions');
path(path,p1) %auxiliary functions

% p1 = genpath('D:\CODES_final\functions');
% path(path,p1) %auxiliary functions
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

YA=nanmean(GRIDS,2);

%Construct different periods
GRIDS_1880_2022 = GRIDS(:,all(~isnan(GRIDS)));
YB=mean(GRIDS_1880_2022,2);

GRIDS1960=GRIDS(81:end,:);
GRIDS_1960_2022 = GRIDS1960(:,all(~isnan(GRIDS1960)));
G1960=mean(GRIDS_1960_2022,2);

GRIDS1920=GRIDS(41:end,:);
GRIDS_1920_2022 = GRIDS1920(:,all(~isnan(GRIDS1920)));
G1920=mean(GRIDS_1920_2022,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ADF UR test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1880-2022
[t,n]=size(GRIDS_1880_2022);
pmax=round(12*(t/100)^0.25);
for i=1:n
   g=GRIDS_1880_2022(:,i);
   [adf,p]=calcula_adf_t_sbic(g,pmax); %ADF with intercep and trend and SBIC
   [h,pValue,stat,cValue,reg]=adftest(g,'model','TS','lags',p,'Alpha', 0.01); %Compute ADF
   H(i)=h;
end
reject_adf_1880_2022=sum(H)/n*100;
MATRIX(1,1)=reject_adf_1880_2022;

%1920-2022
[t,n]=size(GRIDS_1920_2022);
pmax=round(12*(t/100)^0.25);
for i=1:n
   g=GRIDS_1920_2022(:,i);
   [adf,p]=calcula_adf_t_sbic(g,pmax); %ADF with intercep and trend and SBIC
   [h,pValue,stat,cValue,reg]=adftest(g,'model','TS','lags',p,'Alpha', 0.01); %Compute ADF
   H(i)=h;
end
reject_adf_1920_2022=sum(H)/n*100;
MATRIX(2,1)=reject_adf_1920_2022;

%1960-2022
[t,n]=size(GRIDS_1960_2022);
pmax=round(12*(t/100)^0.25);
for i=1:n
   g=GRIDS_1960_2022(:,i);
   [adf,p]=calcula_adf_t_sbic(g,pmax); %ADF with intercep and trend and SBIC
   [h,pValue,stat,cValue,reg]=adftest(g,'model','TS','lags',p,'Alpha', 0.01); %Compute ADF
   H(i)=h;
end
reject_adf_1960_2022=sum(H)/n*100;
MATRIX(3,1)=reject_adf_1960_2022;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       UR-SB test Kim-Perron (2009)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1880-2022
[t,n]=size(GRIDS_1880_2022);
for i=1:n
   g=GRIDS_1880_2022(:,i);
   Dt=breakestimate(g,3); %estimate break, model 3
   khat=MIC(g,0,2,Dt); %select lag with MBIC and model with intercept and slope
   test=UR3k(g,Dt,khat); %Unit root test
   lambda=round(Dt/t,1); sig=1;
     if lambda==0
        lambda=0.1;
     elseif lambda==1
        lambda=0.9;
    end
    cv=TableIVB_Perron_1989(lambda,sig);
    if abs(test)>abs(cv)
        UR(i)=0;
    else
        UR(i)=1;
    end
end
reject_kp_1880_2022=(n-sum(UR))/n*100;
MATRIX(1,2)=reject_kp_1880_2022;

%1920-2022
[t,n]=size(GRIDS_1920_2022);
for i=1:n
   g=GRIDS_1920_2022(:,i);
   Dt=breakestimate(g,3); %estimate break, model 3
   khat=MIC(g,0,2,Dt); %select lag with MBIC and model with intercept and slope
   test=UR3k(g,Dt,khat); %Unit root test
   lambda=round(Dt/t,1); sig=1;
     if lambda==0
        lambda=0.1;
     elseif lambda==1
        lambda=0.9;
    end
    cv=TableIVB_Perron_1989(lambda,sig);
    if abs(test)>abs(cv)
        UR(i)=0;
    else
        UR(i)=1;
    end
end
reject_kp_1920_2022=(n-sum(UR))/n*100;
MATRIX(2,2)=reject_kp_1920_2022;

%1960-2022
[t,n]=size(GRIDS_1960_2022);
for i=1:n
   g=GRIDS_1960_2022(:,i);
   Dt=breakestimate(g,3); %estimate break, model 3
   khat=MIC(g,0,2,Dt); %select lag with MBIC and model with intercept and slope
   test=UR3k(g,Dt,khat); %Unit root test
   lambda=round(Dt/t,1); sig=1;
     if lambda==0
        lambda=0.1;
     elseif lambda==1
        lambda=0.9;
    end
    cv=TableIVB_Perron_1989(lambda,sig);
    if abs(test)>abs(cv)
        UR(i)=0;
    else
        UR(i)=1;
    end
end
reject_kp_1960_2022=(n-sum(UR))/n*100;
MATRIX(3,2)=reject_kp_1960_2022;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SB Perron_Yabu (2009)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model=3;  
criteria=2; %BIC if criteria=2, AIC if criteria=1
eps=0.15;	%chose a value from {0.01, 0.05, 0.10, 0.15, 0.25}

%1880-2022
[t,n]=size(GRIDS_1880_2022);
kmax=fix(12*(t/100)^(1/4)); %maximum number of lag length
for i=1:n
    g=GRIDS_1880_2022(:,i);
    [test,cv,tb]=qfgls_simus(g,kmax,model,criteria,eps);
    if test>cv(3)
        SB(i)=1;
    else
        SB(i)=0;
    end
end

struct_break_1880_2022=sum(SB)/n*100;
MATRIX(1,3)=struct_break_1880_2022;

%1920-2022
[t,n]=size(GRIDS_1920_2022);
kmax=fix(12*(t/100)^(1/4)); %maximum number of lag length
for i=1:n
    g=GRIDS_1920_2022(:,i);
    [test,cv,tb]=qfgls_simus(g,kmax,model,criteria,eps);
    if test>cv(3)
        SB(i)=1;
    else
        SB(i)=0;
    end
end

struct_break_1920_2022=sum(SB)/n*100;
MATRIX(2,3)=struct_break_1920_2022;

%1960-2022
[t,n]=size(GRIDS_1960_2022);
kmax=fix(12*(t/100)^(1/4)); %maximum number of lag length
for i=1:n
    g=GRIDS_1960_2022(:,i);
    [test,cv,tb]=qfgls_simus(g,kmax,model,criteria,eps);
    if test>cv(3)
        SB(i)=1;
    else
        SB(i)=0;
    end
end

struct_break_1960_2022=sum(SB)/n*100;

%Data for Table 2
MATRIX(3,3)=struct_break_1960_2022;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Trend parameters linear trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1880-2022
[t,n]=size(GRIDS_1880_2022);
trend=(1:t)';
X=[ones(t,1),trend];
for i=1:n
    g=GRIDS_1880_2022(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    BETA_1880_2022(i)=b(2);
end

%1920-2022
[t,n]=size(GRIDS_1920_2022);
trend=(1:t)';
X=[ones(t,1),trend];
for i=1:n
    g=GRIDS_1920_2022(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    BETA_1920_2022(i)=b(2);
end

%1960-2022
[t,n]=size(GRIDS_1960_2022);
trend=(1:t)';
X=[ones(t,1),trend];
for i=1:n
    g=GRIDS_1960_2022(:,i);
    [b,t_nw, ~, res, ~, ~, ~, ~]=ols_hac(g,X); %Estimate trend
    BETA_1960_2022(i)=b(2);
end

figure(1)
[fx,x]=ksdensity(BETA_1880_2022)
plot(x,fx,'LineWidth',2);
hold on
[fx,x]=ksdensity(BETA_1920_2022)
plot(x,fx,'LineWidth',2);
hold on
[fx,x]=ksdensity(BETA_1960_2022)
plot(x,fx,'LineWidth',2);
hold off
legend('1880-2022','1920-2022','1960-2022')

set( figure(1),'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Figure_3a'
print -dpng 'Figure_3a'
print -deps 'Figure_3a' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Trend parameters with broken trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1880-2022
model=3;  
criteria=2; %BIC if criteria=2, AIC if criteria=1
eps=0.15;	%chose a value from {0.01, 0.05, 0.10, 0.15, 0.25}

[t,n]=size(GRIDS_1880_2022);
pmax=round(12*(t/100)^0.25);
for i=1:n
    g=GRIDS_1880_2022(:,i);
    [test,cv,tb]=qfgls_simus(g,pmax,model, criteria,eps);
    TB(i)=tb;
end
trend=(1:t)';
for i=1:n
    g=GRIDS_1880_2022(:,i);
    tb=TB(i);
    Du=[zeros(tb,1);ones(t-tb,1)];
    Dt=[zeros(tb,1);trend(tb+1:end)-tb];
    X=[ones(t,1),Du,trend,Dt];
    [b]=ols_hac(g,X); %Estimate model with SB
    GAMMA1(i)=b(3);
    GAMMA2(i)=b(4);
end

figure(2)
[fx,x]=ksdensity(GAMMA1);
plot(x,fx,"LineWidth",2)
hold on
[fx,x]=ksdensity(GAMMA1+GAMMA2);
plot(x,fx,"LineWidth",2)
legend('$\gamma_{1}$','$\gamma_{1}+\gamma_{2}$','Interpreter','latex')

set( figure(2),'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Figure_3b'
print -dpng 'Figure_3b'
print -deps 'Figure_3b' 