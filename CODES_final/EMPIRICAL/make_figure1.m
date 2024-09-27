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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Working with all grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G=NaN(nyears-1,ns);
for i=1:nyears-1
    Y=DATA((i-1)*12+1:i*12,:);
    G(i,:)=nanmean(Y);
end

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


Gsh=NaN(nyears-1,ns2);
for i=1:nyears-1
    Ysh=Temp_sh((i-1)*12+1:i*12,:);
    Gsh(i,:)=nanmean(Ysh);
end


for i=1:nyears-2
    G2=G(i:end,:);
    G2= G2(:,all(~isnan(G2)));
    N(i)=size(G2,2);
end

for i=1:nyears-2
    G2=Gnh(i:end,:);
    G2= G2(:,all(~isnan(G2)));
    Nnh(i)=size(G2,2);
end

for i=1:nyears-2
    G2=Gsh(i:end,:);
    G2= G2(:,all(~isnan(G2)));
    Nsh(i)=size(G2,2);
end

years=years(1:end-2);
figure(1)
plot(years,N,'LineWidth',2)
hold on
plot(years,Nnh,'LineWidth',2)
hold on
plot(years,Nsh,'LineWidth',2)
hold off
xlim([1850 2021]);
legend('Globe','North Hemisphere','South Hemisphere','Location','best')

set( figure(1),'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Figure_1a'
print -dpng 'Figure_1a'
print -deps 'Figure_1a' 

for i=1:nyears-2
    NM(i)=(ns-sum(isnan(G(i,:))))/ns;
end


for i=1:nyears-2
    NMnh(i)=(ns2-sum(isnan(Gnh(i,:))))/ns2;
end

for i=1:nyears-2
    NMsh(i)=(ns2-sum(isnan(Gsh(i,:))))/ns2;
end

figure(2)
plot(years,NM,'LineWidth',2)
hold on
plot(years,NMnh,'LineWidth',2)
hold on
plot(years,NMsh,'LineWidth',2)
hold off
xlim([1850 2021]);
legend('Globe','North Hemisphere','South Hemisphere','Location','best')

set( figure(2),'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Figure_2a'
print -dpng 'Figure_2a'
print -deps 'Figure_2a' 