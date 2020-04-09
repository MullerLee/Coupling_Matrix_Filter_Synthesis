% filename: Cross_5_Proposed_Poles.m
% author:   Zve.L
% date:     3/11/2020  
clear all;
syms u1 v1 u2 v2 u3 v3 u4 v4 u5 v5 w w0
%Parameters
fL=3300; fH=3700; f0=sqrt(fL*fH);
FBW=(fH-fL)/f0;
Pol1=3900; Pol2=4000; Pol3=4200;
w1=(Pol1/f0-f0/Pol1)/FBW;
w2=(Pol2/f0-f0/Pol2)/FBW;
w3=(Pol3/f0-f0/Pol3)/FBW;
w0=sqrt(w^2-1);
% Degree 1
u1=-1/w1+w;
v1=w0*sqrt(1-1/(w1^2));
% Degree 2
u2=w*u1-u1/w2+w0*v1*sqrt(1-1/(w2^2));
v2=w*v1-v1/w2+w0*u1*sqrt(1-1/(w2^2));
% Degree 3
u3=w*u2-u2/w3+w0*v2*sqrt(1-1/(w3^2));
v3=w*v2-v2/w3+w0*u2*sqrt(1-1/(w3^2));
% Degree 4
u4=w*u3+w0*v3;
v4=w*v3+w0*u3;
% Degree 5
u5=w*u4+w0*v4;
v5=w*v4+w0*u4;

u5c=sym2poly(u5);
S11_Zeros=double(solve(u5));
%u5=u5/u5c(1);
%u5c=sym2poly(u5);