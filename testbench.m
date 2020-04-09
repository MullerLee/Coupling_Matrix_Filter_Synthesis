count=1;

RL = 20; s = 1i;
F_val = polyval(sym2poly(F),s);
P_val = polyval(sym2poly(P),s);
E_val = polyval(sym2poly(E),s);
e=P_val/(F_val*sqrt(10^(RL/10)-1));

for x0=-3:0.01:5
    s=x0*1i;
    F_val = polyval(sym2poly(F),s);
    P_val = polyval(sym2poly(P),s);
    E_val = polyval(sym2poly(E),s);
    S21(count)=-10*log(sqrt(1/(1+(e^2)*(F_val/P_val)^2)));
    S11(count)=-10*log(sqrt((e^2)*(F_val/P_val)^2/(1+(e^2)*(F_val/P_val)^2)));
    S112(count)=-10*log(F_val/E_val);
    count=count+1;
    disp(count);
end

xx=-3:0.01:5;
plot(xx,S112,'r','linewidth',2);hold on;
plot(xx,S11,'b','linewidth',2);
%hold off;
grid on;
set(gca,'linewidth',2)
xlabel('LOWPASS PROTOTYPE FREQUENCY (rad/dec)','fontsize',14)
ylabel('RETURN LOSS (dB)','fontsize',14)
legend('S21','S11');