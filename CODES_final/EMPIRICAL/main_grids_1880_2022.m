clear 
%Build means with method A and B with grids data from 1880-2022
%Compute ADF-UR test
p1 = genpath('d:\WORKA\Climate\UR_letter\CODES_final\functions');
path(path,p1) %auxiliary functions
load DATA_GRID
DATA=DATA_GRID.temp;
years=unique(DATA_GRID.years);
nyears=length(years);
ns=size(DATA,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Working with all grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G=NaN(nyears-1,ns);
for i=1:nyears-1
    Y=DATA((i-1)*12+1:i*12,:);
    G(i,:)=nanmean(Y);
end
G=G(31:end,:); %start in 1880, end in 2022
T=size(G,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Reproduce global mean with all grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=nanmean(G,2);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  By hemispheres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAT=DATA_GRID.LAT;
index_nh=LAT>=0;

Temp_nh=DATA(:,index_nh==1);
Temp_sh=DATA(:,index_nh==0);

ns2=ns/2;

Gnh=NaN(nyears-1,ns2);
for i=1:nyears-1
    Ynh=Temp_nh((i-1)*12+1:i*12,:);
    Gnh(i,:)=nanmean(Ynh);
end
Gnh=Gnh(31:end,:); %start in 1880, end in 2022

znh=nanmean(Gnh,2);



Gsh=NaN(nyears-1,ns2);
for i=1:nyears-1
    Ysh=Temp_sh((i-1)*12+1:i*12,:);
    Gsh(i,:)=nanmean(Ysh);
end
Gsh=Gsh(31:end,:); %start in 1880, end in 2022

zsh=nanmean(Gsh,2);

z2=(2*znh+zsh)/3; %compute mean according with CRU

figure(1)
time=years(31:end-1);
plot(time,z2,time,znh,time,zsh,'LineWidth',2)
xlim([1880 2022])
legend('Globe','North Hemisphere', 'South Hemisphere','Location','best')
title('Temperature average (A)')
set( figure(1),'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Figure_grids_mean_all_grids_1880_2022'
print -dpng 'Figure_grids_mean_all_grids_1880_2022'
print -deps 'Figure_grids_mean_all_grids_1880_2022' 

%%%%%%%Test UR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmax=round(12*(T/100)^0.25);
[adf_z,p]=calcula_adf_t_sbic(z2,pmax); %ADF with intercep and trend and SBIC
[h,pValue,stat,cValue,reg]=adftest(z2,'model','TS','lags',p,'Alpha', 0.01);
Pvalue_z=pValue;
[adf_nh,p]=calcula_adf_t_sbic(znh,pmax); %ADF with intercep and trend and SBIC
[h,pValue,stat,cValue,reg]=adftest(znh,'model','TS','lags',p,'Alpha', 0.01);
Pvalue_znh=pValue;
[adf_sh,p]=calcula_adf_t_sbic(zsh,pmax); %ADF with intercep and trend and SBIC
[h,pValue,stat,cValue,reg]=adftest(zsh,'model','TS','lags',p,'Alpha', 0.01);
Pvalue_zsh=pValue;

disp('ADF test for global average with all grids')
adf_z
disp('pvalue of ADF for global average with all grids')
Pvalue_z

disp('ADF test for NH average with all grids')
adf_nh
disp('pvalue of ADF for NH average with all grids')
Pvalue_znh

disp('ADF test for SH average with all grids')
adf_sh
disp('pvalue of ADF for SH average with all grids')
Pvalue_zsh

DATA_GRID.methodA.global_mean=z2;
DATA_GRID.methodA.nh_mean=znh;
DATA_GRID.methodA.sh_mean=zsh;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Working with selected grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load grids_yearly_selected_1880_2022
[t,n]=size(G_selected.data);


trend=(1:t)';
X=[ones(t,1),trend];
for i=1:n
    g=G_selected.data(:,i);
    [b,t_nw, se_nw, res, r2, varres, varBhat, y_hat]=ols_hac(g,X); %Estimate trend
    BETA(i)=b(2);
    [phi, t_nw, se_nw,res, r2, varres, varBhat, y_hat]=arols(res,1,1); %Estimate ar 
    AR(i)=phi(2);
    [adf,p]=calcula_adf_t_sbic(g,pmax); %ADF with intercep and trend and SBIC
    ADF(i)=adf;
    [h,pValue,stat,cValue,reg]=adftest(g,'model','TS','lags',p,'Alpha', 0.01); %Compute ADF
    if pValue<=0.05
        UR(i)=0;
    else
        UR(i)=1;
    end
end
disp('The % of grids with UR is:')
sum(UR)/n


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Analizing mean with stable grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_selected_grids=mean(G_selected.data,2);
[adf,p]=calcula_adf_t_sbic(z_selected_grids,pmax); %ADF with intercep and trend and SBIC
[h,pValue,stat,cValue,reg]=adftest(z_selected_grids,'model','TS','lags',p,'Alpha', 0.01); %Compute ADF

disp('ADF test for the Globe with selected grids')
adf
disp('pvalue of ADF for the Globe with seledted grids')
pValue

index_nh=G_selected.lat>=0;
z_nh_selected_grids=mean(G_selected.data(:,index_nh==1),2);
z_sh_selected_grids=mean(G_selected.data(:,index_nh==0),2);

[adf_nh,p]=calcula_adf_t_sbic(z_nh_selected_grids,pmax); %ADF with intercep and trend and SBIC
[h,pValue_nh,stat,cValue,reg]=adftest(z_nh_selected_grids,'model','TS','lags',p,'Alpha', 0.01); %Compute ADF

disp('ADF test for the NH with selected grids')
adf_nh
disp('pvalue of ADF for the NH with seledted grids')
pValue_nh

[adf_sh,p]=calcula_adf_t_sbic(z_sh_selected_grids,pmax); %ADF with intercep and trend and SBIC
[h,pValue_sh,stat,cValue,reg]=adftest(z_sh_selected_grids,'model','TS','lags',p,'Alpha', 0.01); %Compute ADF

disp('ADF test for the SH with selected grids')
adf_sh
disp('pvalue of ADF for the SH with seledted grids')
pValue_sh

DATA_GRID.methodB.global_mean=z_selected_grids;
DATA_GRID.methodB.nh_mean=z_nh_selected_grids;
DATA_GRID.methodB.sh_mean=z_sh_selected_grids;

figure(2)
time=years(31:end-1);
plot(time,z_selected_grids,time,z_nh_selected_grids,time,z_sh_selected_grids,'LineWidth',2)
xlim([1880 2022])
legend('Globe','North Hemisphere', 'South Hemisphere','Location','best')
title('Temperature average (B)')
set( figure(2),'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Figure_grids_mean_selected_grids_1880_2022'
print -dpng 'Figure_grids_mean_selected_grids_1880_2022'
print -deps 'Figure_grids_mean_selected_grids_1880_2022' 

save('DATA_GRID_1880_2022.mat','DATA_GRID')
