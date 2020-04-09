% filename: Couping_Matrix_Extraction.m
% author:   Zve.L
% date:     3/12/2020
% rev.log   
% Please Run After Cross_5_Proposed_Poles.m

% extracting epsilon
RL = 5; s=1i;
% P=(1-s/w1)*(1-s/w2)*(1-s/w3);
% F=u5c(6)+u5c(5)*s+u5c(4)*s^2+u5c(3)*s^3+u5c(2)*s^4+u5c(1)*s^5;
P=(1-s/(1i*w1))*(1-s/(1i*w2))*(1-s/(1i*w3));
F=(1-s/(1i*S11_Zeros(1)))*(1-s/(1i*S11_Zeros(2)))*(1-s/(1i*S11_Zeros(3)))*(1-s/(1i*S11_Zeros(4)))*(1-s/(1i*S11_Zeros(5)));
F=F*(-(1i)^5*S11_Zeros(1)*S11_Zeros(2)*S11_Zeros(3)*S11_Zeros(4)*S11_Zeros(5)); 
P=P*(-(1i)^3*w1*w2*w3);                      %%% Unity the coefficient
e=P/(F*sqrt(10^(RL/10)-1));

% polyminal operation, transfer to s-plane
syms s S21 E2 E;
P=(1-s/(1i*w1))*(1-s/(1i*w2))*(1-s/(1i*w3));
%F=u5c(6)+u5c(5)*(s/1i)+u5c(4)*(s/1i)^2+u5c(3)*(s/1i)^3+u5c(2)*(s/1i)^4+u5c(1)*(s/1i)^5;
F=(1-s/(1i*S11_Zeros(1)))*(1-s/(1i*S11_Zeros(2)))*(1-s/(1i*S11_Zeros(3)))*(1-s/(1i*S11_Zeros(4)))*(1-s/(1i*S11_Zeros(5)));
E2=((e^2)*(F^2)+(P^2))/(e^2);
F=F*(-(1i)^5*S11_Zeros(1)*S11_Zeros(2)*S11_Zeros(3)*S11_Zeros(4)*S11_Zeros(5)); 
P=P*(-(1i)^3*w1*w2*w3);                      %%% Unity the coefficient
% P=(1-s/(w1))*(1-s/(w2))*(1-s/(w3));
% F=u5c(6)+u5c(5)*(s)+u5c(4)*(s)^2+u5c(3)*(s)^3+u5c(2)*(s)^4+u5c(1)*(s)^5;
% E2=((e^2)*(F^2)+(P^2))/(e^2);
% E=sqrt(E2);

Eroot=double(solve(E2));
Ecof=Eroot(1:5);

E=(1-s/Ecof(1))*(1-s/Ecof(2))*(1-s/Ecof(3))*(1-s/Ecof(4))*(1-s/Ecof(5));
%E=E*(-Ecof(1)*Ecof(2)*Ecof(3)*Ecof(4)*Ecof(5));  %%% Unity the coefficient
Ecf=sym2poly(E);

% Build [y], since N is odd
syms y21 y11 m1 n1;
m1=real(Ecf(1))+imag(Ecf(2))*s+real(Ecf(3))*s^2+imag(Ecf(4))*s^3+real(Ecf(5))*s^4+imag(Ecf(6))*s^5;
n1=imag(Ecf(1))+real(Ecf(2))*s+imag(Ecf(3))*s^2+real(Ecf(4))*s^3+imag(Ecf(5))*s^4+real(Ecf(6))*s^5;
y21 = P/(e*n1);
y22 = m1/n1;

% Matrix Synthesis
 % residue of [y]
 %syms u1 v1 p1
% %  [nu1 , dem1]=numden(y22);
% %  [nu2 , dem2]=numden(y21);
% %  [r1,p1,k1] = residue(sym2poly(nu1) , sym2poly(dem1));
% %  [r2,p2,k2] = residue(sym2poly(nu2) , sym2poly(dem2));  % Note that P1 & p2 are identical and ...
% %                                                         % be eigenvals of M
 % Optimized Way
 [r2,p2,k2] = residue(sym2poly(m1) , sym2poly(n1));     % r2 for y22 --r22k
 [r1,p1,k1] = residue(sym2poly(P) , sym2poly(e*n1));    % r1 for y21 --r21k
 
 % fulfill T matrix 1 and N coloum
 TN = sqrt(r2);
 T1 = r1./sqrt(r2);
 para_n1=sqrt(sum(T1.^2));
 para_n2=sqrt(sum(TN.^2));
 T1_o=T1/para_n1;
 TN_o=TN/para_n2;
 T1_o=T1_o';
 TN_o=TN_o';
 T1_o=real(T1_o/norm(T1_o));
 TN_o=TN_o/norm(TN_o); % Basic Norm to eliminate calculation error
 TN_o=TN_o/norm(TN_o);
 
 % fulfill whole Matrix with Gram-Schmidt Orthonormal
 t2=[1,2,3,5,6];
 t3=[6,2,5,6,8];
 t4=[1,5,6,3,5]; % Three random non-linear related vects
 %Ditermin = det([T1_o,t2,t3,t4,TN_o]); % Check if linear orthonormal
 T2=t2-dot(T1_o,t2)/dot(T1_o,T1_o)*T1_o-dot(TN_o,t2)/dot(TN_o,TN_o)*TN_o; %Pick Any one Method
 %T2=t2-((T1_o'*t2)./(T1_o'*T1_o))*T1_o-((TN_o'*t2)./(TN_o'*TN_o))*TN_o;
 T3=t3-dot(T1_o,t3)/dot(T1_o,T1_o)*T1_o-dot(TN_o,t3)/dot(TN_o,TN_o)*TN_o-dot(T2,t3)/dot(T2,T2)*T2; %Pick Any one Method
 T4=t4-dot(T1_o,t4)/dot(T1_o,T1_o)*T1_o-dot(TN_o,t4)/dot(TN_o,TN_o)*TN_o-dot(T2,t4)/dot(T2,T2)*T2-dot(T3,t4)/dot(T3,T3)*T3; %Pick Any one Method
 T2=T2/norm(T2);
 T3=T3/norm(T3);
 T4=T4/norm(T4);
 T_Full=[T1_o;T2;T3;T4;TN_o];  % The Built T Matrix (With Error)
 
% Build Coupling Matrix
  %pre_regulation
  p1=imag(p1);
M=T_Full*diag(p1)*T_Full';
M=real(M);