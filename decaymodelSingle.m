function I = decaymodelSingle(x,xdata)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
%x(1) const bg; 
%x(2) amplitude c1; 
%x(3) tau lifetime; 
%x(4) temporal shift t0; 
I = x(1)+x(2).*exp(-(xdata-x(4))./x(3));
end