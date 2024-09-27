function cv=TableIVB_Perron_1989 (lambda,sig)

LAMBDA= [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
sig1=[-4.30 -4.39 -4.39 -4.34 -4.32 -4.45 -4.42 -4.33, -4.27];
sig5=[-3.68, -3.77 -3.76 -3.72, -3.76 -3.76, -3.80, -3.75, -3.69];
sig10=[-3.40, -3.47 -3.46, -3.44, -3.46 -3.47, -3.51, -3.46, -3.38];
pos=find(LAMBDA==lambda);
if sig==1
    cv=sig1(pos);
elseif sig==5
    cv=sig5(pos);
elseif sig==10
    cv=sig10(pos);
end
