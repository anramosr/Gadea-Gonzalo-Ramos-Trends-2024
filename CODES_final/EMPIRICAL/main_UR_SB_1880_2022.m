close all; clear

%Test for UR with structural breaks with grid data 1880-2022

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

%%%%%%%%%%%%%%%%%%%%%%%%%%Kim and Perron 2009%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The Globe
Dt=breakestimate(z2,3) %estimate break, model 3
khat=MIC(z2,0,2,Dt); %select lag with MBIC and model with intercept and slope
test=UR3k(z2,Dt,khat) %Unit root test
lambda=round(Dt/T,1); sig=1;
cv=TableIVB_Perron_1989 (lambda,sig)

%The NH
Dt=breakestimate(z_nh,3) %estimate break, model 3
khat=MIC(z_nh,0,2,Dt); %select lag with MBIC and model with intercept and slope
test=UR3k(z_nh,Dt,khat) %Unit root test
lambda=round(Dt/T,1); sig=1;
cv=TableIVB_Perron_1989 (lambda,sig)

%The SH
Dt=breakestimate(z_sh,3) %estimate break, model 3
khat=MIC(z_sh,0,2,Dt); %select lag with MBIC and model with intercept and slope
test=UR3k(z_sh,Dt,khat) %Unit root test
lambda=round(Dt/T,1); sig=1;
cv=TableIVB_Perron_1989 (lambda,sig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Working with selected grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load grids_yearly_selected_1880_2022
[T,n]=size(G_selected.data);
for i=1:n
    g=G_selected.data(:,i);
    Dt=breakestimate(g,3); %estimate break, model 3
    SB(i)=Dt;
    khat=MIC(g,0,2,Dt); %select lag with MBIC and model with intercept and slope
    test=UR3k(g,Dt,khat); %Unit root test
    KP(i)=test;
    lambda=round(Dt/T,1); sig=1;
     if lambda==0
        lambda=0.1;
    elseif lambda==1
        lambda=0.9;
    end
    cv=TableIVB_Perron_1989 (lambda,sig);
    CV(i)=cv;
    if abs(test)>abs(cv)
        UR(i)=0;
    else
        UR(i)=1;
    end
end


disp('The % of series with UR is:')
sum(UR)/n

z=DATA_GRID.methodB.global_mean;
z_nh=DATA_GRID.methodB.nh_mean;
z_sh=DATA_GRID.methodB.sh_mean;

%The Globe
Dt=breakestimate(z,3) %estimate break, model 3
khat=MIC(z,0,2,Dt); %select lag with MBIC and model with intercept and slope
test=UR3k(z,Dt,khat) %Unit root test
lambda=round(Dt/T,1); sig=1;
cv=TableIVB_Perron_1989 (lambda,sig)

%The NH
Dt=breakestimate(z_nh,3) %estimate break, model 3
khat=MIC(z_nh,0,2,Dt); %select lag with MBIC and model with intercept and slope
test=UR3k(z_nh,Dt,khat) %Unit root test
lambda=round(Dt/T,1); sig=1;
cv=TableIVB_Perron_1989 (lambda,sig)

%The SH
Dt=breakestimate(z_sh,3) %estimate break, model 3
khat=MIC(z_sh,0,2,Dt); %select lag with MBIC and model with intercept and slope
test=UR3k(z_sh,Dt,khat) %Unit root test
lambda=round(Dt/T,1); sig=1;
cv=TableIVB_Perron_1989 (lambda,sig)