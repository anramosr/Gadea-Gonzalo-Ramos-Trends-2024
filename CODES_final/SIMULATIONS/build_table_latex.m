clear
load RDOS_simus_v2_nrep1000
%The structure correspons with 4 methods of agreggation, 2 DGD, 2 ways to compute the mean

%for Table 3, DGP 1, mean A
mat_adf_1_1=mean(squeeze(RDOS.ADF.reject(:,1,1,:)),2)*100; %ADF mean A, linear model
mat_adf_1_1=100*ones(4,1)-mat_adf_1_1; %convert to no-rejections
mat_kp_1_1=mean(squeeze(RDOS.KP.reject(:,1,1,:)),2)*100; %ADF mean A, linear model
mat_kp_1_1=100*ones(4,1)-mat_kp_1_1; %convert to no-rejections
mat_py_1_1=mean(squeeze(RDOS.PY.reject(:,1,1,:)),2)*100; %ADF mean A, linear model

%for Table 3, DGP 1, mean B
mat_adf_1_2=mean(squeeze(RDOS.ADF.reject(:,1,2,:)),2)*100; %ADF mean B, linear model
mat_adf_1_2=100*ones(4,1)-mat_adf_1_2; %convert to no-rejections
mat_kp_1_2=mean(squeeze(RDOS.KP.reject(:,1,2,:)),2)*100; %ADF mean B, linear model
mat_kp_1_2=100*ones(4,1)-mat_kp_1_2; %convert to no-rejections
mat_py_1_2=mean(squeeze(RDOS.PY.reject(:,1,2,:)),2)*100; %ADF mean B, linear model

Matrix_DPG1=[mat_adf_1_1,mat_kp_1_1,mat_py_1_1,mat_adf_1_2,mat_kp_1_2,mat_py_1_2];

%for Table 4, DGP 2, mean A
mat_adf_2_1=mean(squeeze(RDOS.ADF.reject(:,2,1,:)),2)*100; %ADF mean A, linear model
mat_adf_2_1=100*ones(4,1)-mat_adf_2_1; %convert to no-rejections
mat_kp_2_1=mean(squeeze(RDOS.KP.reject(:,2,1,:)),2)*100; %ADF mean A, linear model
mat_kp_2_1=100*ones(4,1)-mat_kp_2_1; %convert to no-rejections
mat_py_2_1=mean(squeeze(RDOS.PY.reject(:,2,1,:)),2)*100; %ADF mean A, linear model

%for Table 4, DGP 2, mean B
mat_adf_2_2=mean(squeeze(RDOS.ADF.reject(:,2,2,:)),2)*100; %ADF mean B, linear model
mat_adf_2_2=100*ones(4,1)-mat_adf_2_2; %convert to no-rejections
mat_kp_2_2=mean(squeeze(RDOS.KP.reject(:,2,2,:)),2)*100; %ADF mean B, linear model
mat_kp_2_2=100*ones(4,1)-mat_kp_2_2; %convert to no-rejections
mat_py_2_2=mean(squeeze(RDOS.PY.reject(:,2,2,:)),2)*100; %ADF mean B, linear model

Matrix_DPG2=[mat_adf_2_1,mat_kp_2_1,mat_py_2_1,mat_adf_2_2,mat_kp_2_2,mat_py_2_2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Replicate Table 3 of the letter
tableCaption='Table 3: Proportion of non-rejections in the unit root tests and structural break detection analysis on simulated averages (linear-trend model))';
tableLabel='tab-simus1';
file=strcat('Table_simus1_referee3.tex');
fid=fopen(file,'w');
tt=[' ' '\\begin{table}[h!]\\caption{',tableCaption,'}\\label{',tableLabel,'}\\begin{center}\\scalebox{0.7}{\\begin{tabular}{l|cccccc} \\hline \\hline \n'];
fprintf(fid,tt);


tt=['','& \\multicolumn{3}{c}{Method A}' ,'& \\multicolumn{3}{c}{Method B}', '\\\\  \n'];
fprintf(fid,tt);

tt=['Alternative','& ADF', '& KP', '& PY-SB','& ADF', '& KP','& PY-SB''\\\\  \n'];
fprintf(fid,tt);


rowlabels={'1','1*','2', '2*'};

for i=1:size(Matrix_DPG1,1)
    t1=[];
    for j=1:size(Matrix_DPG1,2)
    t1=[t1,'&',num2str(Matrix_DPG1(i,j),'%5.2f') ];  
    end
    tt=[rowlabels{i},t1,'\\\\ \n'];
    fprintf(fid,tt);
end

tt='\\hline \\hline \\end{tabular}}\\end{center}\\end{table}\n';
fprintf(fid,tt);
st=fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Replicate Table 4 of the letter
tableCaption='Proportion of non-rejections in the unit root tests and structural break detection analysis on simulated averages (model with broken-trend)';
tableLabel='tab-simus2';
file=strcat('Table_simus2_referee3.tex');
fid=fopen(file,'w');
tt=[' ' '\\begin{table}[h!]\\caption{',tableCaption,'}\\label{',tableLabel,'}\\begin{center}\\scalebox{0.7}{\\begin{tabular}{l|cccccc} \\hline \\hline \n'];
fprintf(fid,tt);


tt=['','& \\multicolumn{2}{c}{Method A}' ,'& \\multicolumn{2}{c}{Method B}', '\\\\  \n'];
fprintf(fid,tt);

tt=['Alternative','& ADF', '& KP', '& PY-SB','& ADF', '& KP','& PY-SB''\\\\  \n'];
fprintf(fid,tt);


rowlabels={'1','1*','2', '2*'};

for i=1:size(Matrix_DPG2,1)
    t1=[];
    for j=1:size(Matrix_DPG2,2)
    t1=[t1,'&',num2str(Matrix_DPG2(i,j),'%5.2f') ];  
    end
    tt=[rowlabels{i},t1,'\\\\ \n'];
    fprintf(fid,tt);
end

tt='\\hline \\hline \\end{tabular}}\\end{center}\\end{table}\n';
fprintf(fid,tt);
st=fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%Individual series%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rejections in individual series

mean(mean(RDOS.Y1.adf,2))*100
mean(mean(RDOS.Y2.adf,2))*100
mean(mean(RDOS.Y3.adf,2))*100
mean(mean(RDOS.Y4.adf,2))*100

mean(mean(RDOS.Y1.kp,2))*100
mean(mean(RDOS.Y2.kp,2))*100
mean(mean(RDOS.Y3.kp,2))*100
mean(mean(RDOS.Y4.kp,2))*100

%structural breaks in individual series
mean(mean(RDOS.Y1.py,2))*100
mean(mean(RDOS.Y2.py,2))*100
mean(mean(RDOS.Y3.py,2))*100
mean(mean(RDOS.Y4.py,2))*100