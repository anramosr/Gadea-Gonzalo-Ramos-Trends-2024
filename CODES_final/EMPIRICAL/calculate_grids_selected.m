clear 
%Calculate quantiles with monthly data from 1880-2022 or 1960-2022
load DATA_GRID
DATA=DATA_GRID.temp;
LAT=DATA_GRID.LAT
LONG=DATA_GRID.Long;
years=unique(DATA_GRID.years);
nyears=length(years);
ns=size(DATA,2);

B=NaN(nyears-1,ns);
for i=1:nyears-1
    A=DATA((i-1)*12+1:i*12,:);
    B(i,:)=nanmean(A);
end
G=B(31:end,:); %start in 1880, end in 2022
%G=B(111:end,:); %start in 1960, end in 2022
T=size(G,1);

num_g=0;
A=[];
Lat_selected=[];
Long_selected=[];
for i=2:size(G,2)
    i
    g=G(:,i);
    if array_NaN(g)==1
    num_g=num_g+1;
    A=[A,g];
    Lat_selected=[Lat_selected;LAT(i)];
    Long_selected=[Long_selected;LONG(i)];
    end
end
G_selected.data=A;
G_selected.lat=Lat_selected;
G_selected.long=Long_selected;
G_selected.years=years(31:end-1);  %start in 1880, end in 2022
%G_selected.years=years(111:end-1);  %start in 1960, end in 2022


save('grids_yearly_selected_1880_2022.mat','G_selected')
%save('grids_yearly_selected_1960_2022.mat','G_selected')