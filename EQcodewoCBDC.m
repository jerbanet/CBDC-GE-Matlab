clear;
clc;

beta=0.97;
yl=1;
yh=2;
N=5;
gammad=0.06;
gammae=0.0059;
alpha=0.66;
A=1.45;
Rc=0.98;
Re=1.15;
Rr=1.02;
Ksi=0.06;


y=[yl:0.01:yh]; 



syms Rdh
assume(Rdh > 1)
solve(1/(Rdh-Rc)*(log(Rdh)/beta-log(Rc)/beta+gammad)==yh,Rdh);
dRdh=double(ans)
syms Rdl
assume(Rdl > 1)
solve(1/(Rdl-Rc)*(log(Rdl)/beta-log(Rc)/beta+gammad)==yl,Rdl);
dRdl=double(ans)

Rd=[1:.001:2];
yd=1./(Rd-Rc).*(log(Rd)/beta-log(Rc)/beta+gammad);
figure
plot(Rd,yd)
hold on
yline(yl,'--')
yline(yh,'--')
legend('Deposits to fiat indifference income level')
xlabel('R^d')
ylabel('y*(R^d)')
ytxt = {'y^h';char(repmat(' ',6,1))};
hold off

invbeta=1/beta;
D=.5*yh^2-yh*invbeta./Rd-.5*yd.*yd+invbeta./Rd.*yd;

Dact=D;
    for i=1:length(Rd)
        if Rd(i)<dRdh(1)
            Dact(i)=0;
        end
        if Rd(i)>dRdl(1)
            Dact(i)=.5*(yh^2-yl^2)-(yh-yl)*1/(beta*Rd(i));
        end
        
    end

    StepX=0.001;
Minx=min(Dact)+StepX;
Maxx=max(Dact);
Xarray=[Minx:StepX:Maxx];
Yarray=[Minx:StepX:Maxx;Minx:StepX:Maxx;Minx:StepX:Maxx;Minx:StepX:Maxx];
Eqarray=[Minx:StepX:Maxx;Minx:StepX:Maxx;Minx:StepX:Maxx;Minx:StepX:Maxx];
z=1;
for i=1:length(Xarray)
    valuesearch=Xarray(i);
    while Dact(z)<valuesearch
    z=z+1;
    end
    Yarray(1,i)=Dact(z);
    Yarray(2,i)=Dact(z-1);
    Yarray(3,i)=Rd(z);
    Yarray(4,i)=Rd(z-1);
    Yarray(5,i)=Rd(z-1)+(valuesearch-Dact(z-1))*(Rd(z)-Rd(z-1))/(Dact(z)-Dact(z-1));
    Yarray(6,i)=(Rd(z)-Rd(z-1))/(Dact(z)-Dact(z-1));
    
end

for i=1:length(Xarray)
   
Eqarray(1,i)=Yarray(6,i)*Xarray(i)/N;
Eqarray(2,i)=Yarray(6,i)*Xarray(i)/N+Yarray(5,i);
Eqarray(3,i)=A*alpha*((1-Ksi)*Xarray(i))^(alpha-1);
Eqarray(4,i)=(1-Ksi)*Eqarray(3,i)+Ksi*Rr;

    if Eqarray(4,i)-Eqarray(2,i)>0
    resultplace=i;
    end
end;

figure
plot(Xarray,Yarray(5,:),Xarray,Eqarray(2,:),Xarray,Eqarray(4,:))
ExpArray=[Xarray' Yarray(5,:)' Eqarray(2,:)' Eqarray(4,:)'];
xlabel('R^d')
ylabel('Nd')
legend('R^d(Nd)','dR^d(Nd)d+R^d(Nd)','Loan Supply')


resultplace
Ndstar=Xarray(resultplace)
Rdstar=Yarray(5,resultplace)
ydstar=1/(Rdstar-Rc)*(log(Rdstar)/beta-log(Rc)/beta+gammad)

gammae=0.06;
syms Remin
assume(Remin > 1)
solve(1/(Rdstar-Remin)*(log(Rdstar)/beta-log(Remin)/beta+gammad-gammae)==ydstar,Remin);
dRemin=double(ans)
yestar=1/(Re-Rc)*(log(Re)/beta-log(Rc)/beta+gammae)

gammaei=[0.01:.001:0.06];
dRemini=[0.01:.001:0.06];
for i=1:length(gammaei)
    gammae=gammaei(i);
syms Remin
assume(Remin > 1)
solve(1/(Rdstar-Remin)*(log(Rdstar)/beta-log(Remin)/beta+gammad-gammae)==ydstar,Remin);
dRemin=double(ans);
   dRemini(i)=dRemin;
end;
figure
plot(gammaei,dRemini)

