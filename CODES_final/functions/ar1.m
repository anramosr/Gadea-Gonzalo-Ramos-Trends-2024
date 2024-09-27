function [y]=ar1(ro,ut);
T=length(ut);
y(1)=ut(1);
for i=2:T
    y(i)=ro*y(i-1)+ut(i);
end
y=y';