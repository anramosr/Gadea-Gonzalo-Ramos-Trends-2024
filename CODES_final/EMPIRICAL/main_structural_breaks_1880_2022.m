close all; clear
%Test for structural breaks with grid data 1880-2022
p1 = genpath('d:\WORKA\Climate\UR_letter\CODES_final\functions');
path(path,p1) %auxiliary functions

load DATA_GRID_1880_2022
DATA=DATA_GRID.temp;
years=unique(DATA_GRID.years);
years=years(31:end-1,1);
nyears=length(years);
ns=size(DATA,2);
z2=DATA_GRID.methodA.global_mean;
z_nh=DATA_GRID.methodA.nh_mean;
z_sh=DATA_GRID.methodA.sh_mean;
T=length(z2);


%Choose a model from the following:
%Model=1: Y{t}=a0+a1*DU+b0*t+e{t} where DU=1(t>TB)
%Model=2: Y{t}=a0+b0*t+b1*DT+e{t} where DT=1(t>TB)*(t-TB)
%Model=3: Y{t}=a0+a1*DU+b0*t+b1*DT+e{t} where DU=1(t>TB) and DT=1(t>TB)*(t-TB)
model=3;  
criteria=2; %BIC if criteria=2, AIC if criteria=1

eps=0.15;	%chose a value from {0.01, 0.05, 0.10, 0.15, 0.25}

% T size of the sample

kmax=fix(12*(T/100)^(1/4)); %maximum number of lag length

%***************************************
%Call the main function
%***************************************
%Global mean
qfgls(z2,kmax,model, criteria,eps)
%NH mean
qfgls(z_nh,kmax,model, criteria,eps)
%SH mean
qfgls(z_sh,kmax,model, criteria,eps)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Working with selected grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load grids_yearly_selected_1880_2022
[t,n]=size(G_selected.data);

%G_mean=mean(G_selected.data,2);
G_mean=DATA_GRID.methodB.global_mean;
G_mean_nh=DATA_GRID.methodB.nh_mean;
G_mean_sh=DATA_GRID.methodB.sh_mean;

qfgls(G_mean,kmax,model, criteria,eps)
qfgls(G_mean_nh,kmax,model, criteria,eps)
qfgls(G_mean_sh,kmax,model, criteria,eps)

for i=1:n
    [test,cv,tb]=qfgls_simus(G_selected.data(:,i),kmax,model, criteria,eps);
    PY(i)=test;
    CV(i)=cv(2);
    TB(i)=tb;
    if PY(i)>CV(i)
        SB(i)=1;
    else
        SB(i)=0;
    end
end
disp('The number of series with SB is:')
sum(SB)/n

