clear;
clc;

beta=0.97;
yl=1;
yh=2;
N=5;
gammad=0.06;
gammae=0.02;
alpha=0.66;
A=1.45;
Re=1.02;
Rc=0.98;
Rr=1.02;
Ksi=0.06;


y=[yl:0.01:yh]; 
Rd=[1:.001:2];

yde=1./(Rd-Re).*(log(Rd)/beta-log(Re)/beta+gammad-gammae);
ydc=1./(Rd-Rc).*(log(Rd)/beta-log(Rc)/beta+gammad);
yd=max(yde,ydc);
ye=1/(Re-Rc)*(log(Re)/beta-log(Rc)/beta+gammae);



figure
plot(Rd,ydc,Rd,yde)
hold on
yline(yl,'--','yl')
yline(yh,'--','yh')
yline(ye,'-r')

legend('y^d_c(R^d)','y^d_e(R^d)')
xlabel('R^d')
ylabel('y(R^d)')
ytxt = {'y^h';char(repmat(' ',6,1))};
hold off

Dactstore=createDact(Rd,Re,yd,beta,gammad,gammae,yh,yl,Rc)';


figure
plot(Rd,yd)
gammaei=[0:0.01:gammad];
for gammae=0.01:0.01:gammad
yde=1./(Rd-Re).*(log(Rd)/beta-log(Re)/beta+gammad-gammae);
yd=max(yde,ydc);
Dact=createDact(Rd,Re,yd,beta,gammad,gammae,yh,yl,Rc)';
Dactstore=[Dactstore, Dact];
end

Ndstari=[0:0.01:gammad];
Rdstari=[0:0.01:gammad];
ydstarei=[0:0.01:gammad];
ydstarci=[0:0.01:gammad];
yestari=[0:0.01:gammad];


%Ks=(1-Ksi).*Dact;
%Kd=(1/(alpha*(1-Ksi)).*Rd-1/(alpha*(1-Ksi)).*Ksi*Rr).^(1/(alpha-1));
%figure
%plot(Rd,Ks,Rd,Kd)
%legend('Loan Supply','Loan Demand')

save RdDactfile.mat

