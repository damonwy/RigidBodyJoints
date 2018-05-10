plot(T,K,'-',T,V(1)-V,'-',T,K+V(1)-V,'-');
hold on
%plot(T,Kl,'--',T,Vl-Vl(1),'--',T,Kl+Vl-Vl(1),'--');
%grid
hold off
title('(d) Energy','FontSize',12)
xlabel('Time (s)','FontSize',12)
ylabel('Energy (J)','FontSize',12);
legend('kin.','pot.','tot.');
