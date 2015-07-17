%%% plot the tokance
%%% tqdata contains tl and th, the tokance array at two differential bias arrays
%%% Need to input Vdv, dVd, and sita
%%% plot 3 components of the tokance.
load tqdata

Vdv=-0.5:0.1:0.5;
dVd=0.01;
sita=0.05;
dtq=thM-tlM;
tokance=dtq./dVd./sin(sita);  % in S/m^2
figure()
plot(Vdv, tokance(:,1),'r--','linewidth',[2]); hold on;
plot(Vdv, tokance(:,2),'b--','linewidth',[2]); hold on;
%plot(Vdv, tokance(:,3),'g--','linewidth',[2]); hold on;
legend('in p','p p')
set(gca,'fontsize',[20],'linewidth',[2]);
xlabel('V_D [V]');
ylabel('T [A/m^2]')
grid on
print -dtiff temp



