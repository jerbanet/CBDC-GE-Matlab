clear;
clc;
load RdDactfile.mat
for step=1:length(Ndstari)
    


    StepX=0.001;
    Dact=Dactstore(:,step);
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

    figure
    plot(Xarray,Yarray(5,:))
    plot(Xarray,Yarray(6,:))



    for i=1:length(Xarray)
   
    Eqarray(1,i)=Yarray(6,i)*Xarray(i)/N;
    Eqarray(2,i)=Yarray(6,i)*Xarray(i)/N+Yarray(5,i);
    Eqarray(3,i)=A*alpha*((1-Ksi)*Xarray(i))^(alpha-1);
    Eqarray(4,i)=(1-Ksi)*Eqarray(3,i)+Ksi*Rr;

        if Eqarray(4,i)-Eqarray(2,i)>0
        resultplace=i;
        end
    end;

%end
    plot(Xarray,Yarray(5,:),Xarray,Eqarray(2,:),Xarray,Eqarray(4,:))

    ExpArray=[Xarray' Yarray(5,:)' Eqarray(2,:)' Eqarray(4,:)'];

    resultplace
    Ndstari(step)=Xarray(resultplace)
    Rdstar=Yarray(5,resultplace);
    Rdstari(step)=Yarray(5,resultplace)
    gammae=gammaei(step);
    ydstarei(step)=1/(Rdstar-Re)*(log(Rdstar)/beta-log(Re)/beta+gammad-gammae)
    ydstarci(step)=1/(Rdstar-Rc)*(log(Rdstar)/beta-log(Rc)/beta+gammad)

%gammae=0.01;
%syms Remin
%solve(1/(Rdstar-Remin)*(log(Rdstar)/beta-log(Remin)/beta+gammad-gammae)==ydstar,Remin)
%dRemin=double(ans)
%Rc=0.98;
    yestari(step)=1/(Re-Rc)*(log(Re)/beta-log(Rc)/beta+gammae)
end

save RdDactfile.mat

%filename = 'testdata.xlsx';
%writematrix(ExpArray,filename,'Sheet',1)
