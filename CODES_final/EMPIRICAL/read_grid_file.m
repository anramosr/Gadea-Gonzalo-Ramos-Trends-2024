clear
ncdisp('CRUTEM.5.0.1.0.anomalies.nc')
time  = ncread('CRUTEM.5.0.1.0.anomalies.nc','time');
dates=datevec(time); years=dates(:,1)+1850; months=dates(:,2);
temp  = ncread('CRUTEM.5.0.1.0.anomalies.nc','tas');
latitude=ncread('CRUTEM.5.0.1.0.anomalies.nc','latitude'); 
longitude=ncread('CRUTEM.5.0.1.0.anomalies.nc','longitude'); 
[n,m,T]=size(temp);
G=NaN(T,n*m);
for i=1:T
    G(i,:)=reshape(temp(:,:,i),n*m,1);
end
Long=[];
Lat=[];
for i=1:m
    Lat=[Lat;repmat(latitude(i),n,1)];
end
for i=1:n
    Long=[Long;repmat(longitude(i),m,1)];
end

DATA_GRID.temp=G;
DATA_GRID.years=years;
DATA_GRID.months=months;
DATA_GRID.LAT=Lat;
DATA_GRID.Long=Long;

save('DATA_GRID.mat','DATA_GRID')