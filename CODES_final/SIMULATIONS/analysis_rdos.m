clear 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Four methods of agregation:
%    1 Groups of series enter from time to time randomly
%    1* Series with steeper trends enter later on the record
%    2 Observation process follow a transition probability
%    2* Observation process follow a transition probability
% Two DGP,   Case A: Linear trends, Table 3 paper
%            Case B: Series with Structural Breaks, Table 4 paper
% Two methods to construct the mean:
%            Like CRU, all series
%            With selected grids or stations
% The structure of RDOS. is:
%           ADF.reject (4,2,2,numrep) 
%           ADF.pvalue (4,2,2,numrep) 
% In both cases correspond with:
%           Methods of agregation
%           DGP
%           Type of mean
%           Number of replications

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load RDOS_simus_1000_v1
%% Analysis results UR
% Method agregation 1, DGP A
disp('The percentage of rejections of UR hypothesis with method 1, DGP A for mean_1 is:')
mean(RDOS.ADF.reject(1,1,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 1, DGP A for mean_2 is:')
mean(RDOS.ADF.reject(1,1,2,:))*100

% Method agregation 1, DGP B
disp('The percentage of rejections of UR hypothesis with method 1, DGP B for mean_1 is:')
mean(RDOS.ADF.reject(1,2,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 1, DGP B for mean_2 is:')
mean(RDOS.ADF.reject(1,2,2,:))*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method agregation 1*, DGP A
disp('The percentage of rejections of UR hypothesis with method 1*, DGP A for mean_1 is:')
mean(RDOS.ADF.reject(2,1,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 1*, DGP A for mean_2 is:')
mean(RDOS.ADF.reject(2,1,2,:))*100

% Method agregation 2, DGP B
disp('The percentage of rejections of UR hypothesis with method 1*, DGP B for mean_1 is:')
mean(RDOS.ADF.reject(2,2,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 1*, DGP B for mean_2 is:')
mean(RDOS.ADF.reject(2,2,2,:))*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method agregation 2, DGP A
disp('The percentage of rejections of UR hypothesis with method 2, DGP A for mean_1 is:')
mean(RDOS.ADF.reject(3,1,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 2, DGP A for mean_2 is:')
mean(RDOS.ADF.reject(3,1,2,:))*100

% Method agregation 2, DGP B
disp('The percentage of rejections of UR hypothesis with method 2, DGP B for mean_1 is:')
mean(RDOS.ADF.reject(3,2,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 2, DGP B for mean_2 is:')
mean(RDOS.ADF.reject(3,2,2,:))*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method agregation 2*, DGP A
disp('The percentage of rejections of UR hypothesis with method 2*, DGP A for mean_1 is:')
mean(RDOS.ADF.reject(4,1,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 2*, DGP A for mean_2 is:')
mean(RDOS.ADF.reject(4,1,2,:))*100

% Method agregation 2,* DGP B
disp('The percentage of rejections of UR hypothesis with method 2*, DGP B for mean_1 is:')
mean(RDOS.ADF.reject(4,2,1,:))*100

disp('The percentage of rejections of UR hypothesis with method 2*, DGP B for mean_2 is:')
mean(RDOS.ADF.reject(4,2,2,:))*100


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Analysis structural breaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Analysis structural breaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method agregation 1, DGP A
disp('The percentage of structural breaks with method 1, DGP A for mean_1 is:')
mean(RDOS.PY.reject(1,1,1,:))*100

disp('The percentage of structural breaks with method 1, DGP A for mean_2 is:')
mean(RDOS.PY.reject(1,1,2,:))*100

% Method agregation 1, DGP B
disp('The percentage of structural breaks with method 1, DGP B for mean_1 is:')
mean(RDOS.PY.reject(1,2,1,:))*100

disp('The percentage of structural breaks with method 1, DGP B for mean_2 is:')
mean(RDOS.PY.reject(1,2,2,:))*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method agregation 1*, DGP A
disp('The percentage of structural breaks with method 1*, DGP A for mean_1 is:')
mean(RDOS.PY.reject(2,1,1,:))*100

disp('The percentage of structural breaks with method 1*, DGP A for mean_2 is:')
mean(RDOS.PY.reject(2,1,2,:))*100

% Method agregation 2, DGP B
disp('The percentage of structural breaks with method 1*, DGP B for mean_1 is:')
mean(RDOS.PY.reject(2,2,1,:))*100

disp('The percentage of structural breaks with method 1*, DGP B for mean_2 is:')
mean(RDOS.PY.reject(2,2,2,:))*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method agregation 2, DGP A
disp('The percentage of structural breaks with method 2, DGP A for mean_1 is:')
mean(RDOS.PY.reject(3,1,1,:))*100

disp('The percentage of structural breaks with method 2, DGP A for mean_2 is:')
mean(RDOS.PY.reject(3,1,2,:))*100

% Method agregation 2, DGP B
disp('The percentage of structural breaks with method 2, DGP B for mean_1 is:')
mean(RDOS.PY.reject(3,2,1,:))*100

disp('The percentage of structural breaks with method 2, DGP B for mean_2 is:')
mean(RDOS.PY.reject(3,2,2,:))*100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method agregation 2*, DGP A
disp('The percentage of structural breaks with method 2*, DGP A for mean_1 is:')
mean(RDOS.PY.reject(4,1,1,:))*100

disp('The percentage of structural breaks with method 2*, DGP A for mean_2 is:')
mean(RDOS.PY.reject(4,1,2,:))*100

% Method agregation 2,* DGP B
disp('The percentage of structural breaks with method 2*, DGP B for mean_1 is:')
mean(RDOS.PY.reject(4,2,1,:))*100

disp('The percentage of structural breaks with method 2*, DGP B for mean_2 is:')
mean(RDOS.PY.reject(4,2,2,:))*100

