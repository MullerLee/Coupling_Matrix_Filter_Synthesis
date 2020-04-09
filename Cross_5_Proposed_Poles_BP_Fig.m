% filename: Cross_5_Proposed_Poles_BP_Fig.m
% author:   Zve.L
% date:     3/11/2020 
%Please Run After Cross_5_Proposed_Poles.m, and clean S21, S11
RL = 20; x = 1i;
P=(1-x/(1i*w1))*(1-x/(1i*w2))*(1-x/(1i*w3));
F=(1-x/(1i*S11_Zeros(1)))*(1-x/(1i*S11_Zeros(2)))*(1-x/(1i*S11_Zeros(3)))*(1-x/(1i*S11_Zeros(4)))*(1-x/(1i*S11_Zeros(5)));
F=F*(-(1i)^5*S11_Zeros(1)*S11_Zeros(2)*S11_Zeros(3)*S11_Zeros(4)*S11_Zeros(5)); 
P=P*(-(1i)^3*w1*w2*w3);                      %%% Unity the coefficient
e=P/(F*sqrt(10^(RL/10)-1));
count=1; syms z;
for x0=-3:0.05:5
%     F=u5c(6)+u5c(5)*x+u5c(4)*x^2+u5c(3)*x^3+u5c(2)*x^4+u5c(1)*x^5;
%     P=(1-x/w1)*(1-x/w2)*(1-x/w3);
    x=x0*1i;
    P=(1-x/(1i*w1))*(1-x/(1i*w2))*(1-x/(1i*w3));
    F=(1-x/(1i*S11_Zeros(1)))*(1-x/(1i*S11_Zeros(2)))*(1-x/(1i*S11_Zeros(3)))*(1-x/(1i*S11_Zeros(4)))*(1-x/(1i*S11_Zeros(5)));
    F=F*(-(1i)^5*S11_Zeros(1)*S11_Zeros(2)*S11_Zeros(3)*S11_Zeros(4)*S11_Zeros(5)); 
    P=P*(-(1i)^3*w1*w2*w3);                      %%% Unity the coefficient
    C=F/P;
    buf=1/(1+(e^2)*(C^2));
    S21(count)=10*log(real(sqrt(buf)));
    S11(count)=10*log(real(sqrt(1-buf)));
    sol=solve(z^2-FBW*x0*z*f0-f0^2);
    xx(count)= sol(2);
    count=count+1;
    disp(x);
end

plot(xx,S21,'r','linewidth',2);hold on;
plot(xx,S11,'b','linewidth',2);
%hold off;
grid on
set(gca,'linewidth',2)
xlabel('Freq (Hz)','fontsize',14)
ylabel('Magnitude (dB)','fontsize',14)
legend('S21','S11');