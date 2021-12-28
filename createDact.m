function [Dactres]= createDact(Rd,Re,yd,beta,gammad,gammae,yh,yl,Rc)
if ~exist('Rc','var')
     % third parameter does not exist, so default it to something
      Rc = 0.1;
 end

syms Rdh
assume(Rdh > 1)
solve(1/(Rdh-Re)*(log(Rdh)/beta-log(Re)/beta+gammad-gammae)==yh,Rdh);
dRdh=double(ans);

syms Rdl
assume(Rdl > 1)
solve(1/(Rdl-Re)*(log(Rdl)/beta-log(Re)/beta+gammad-gammae)==yl,Rdl);
dRdle=double(ans);

syms Rdl
assume(Rdl > 1)
solve(1/(Rdl-Rc)*(log(Rdl)/beta-log(Rc)/beta+gammad)==yl,Rdl);
dRdlc=double(ans);
dRdl=max(dRdle,dRdlc);

invbeta=1/beta;
D=.5*yh^2-yh*invbeta./Rd-.5*yd.*yd+invbeta./Rd.*yd;
Dactres=D;
    for i=1:length(Rd)
        if Rd(i)<dRdh
            Dactres(i)=0;
        end
        if Rd(i)>dRdl(1)
            Dactres(i)=.5*(yh^2-yl^2)-(yh-yl)*1/(beta*Rd(i));
        end
        
    end
end