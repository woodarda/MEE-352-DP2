function [CF,TF] = CrossFlow(UA,CF,TF,sys)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
CF.vf=CF.vf(2); CF.P=CF.P(2); CF.hfg=CF.hfg(2); CF.cpf=CF.cpf(2);
CF.muf=CF.muf(2);CF.kf=CF.kf(2);CF.Prf=CF.Prf(2);CF.Rho=CF.Rho(2);
sys=renamevars(sys,"L_g","L");
TF.mdot=TF.mdot(2);
CF.vel=CF.mdot/(CF.Rho*sys.A);
CF.velmax=sys.ST/(sys.ST-sys.D)*CF.vel;
CF.REd=CF.Rho*CF.velmax*sys.D/CF.muf;

if sys.Columns<20
C2=[1 .7 .64; 2 .8 .76;3 .86 .84; 4 .9 .89; 5 .92 .92; 7 .95 .95;...
    10 .97 .97; 13 .98 .98; 16 .99 .99; 17 1 1; 17 1 1;18 1 1;19 1 1];
I=find(C2(:,1)==sys.Columns);
    if sys.Aligned==true
        n=2;
        C2=C2(I,n);
    end
    if sys.Aligned==false
        n=3;
        C2=C2(I,n);
    end
if isempty(I)==true
    C2=[1 .7 .64; 2 .8 .76;3 .86 .84; 4 .9 .89; 5 .92 .92; 7 .95 .95;...
    10 .97 .97; 13 .98 .98; 16 .99 .99; 17 1 1; 17 1 1;18 1 1;19 1 1];
    I=find(C2(:,1)>sys.Columns,1);
    if sys.Aligned ==true
        n=2;
        C2=(sys.Columns-C2(I,1))/(C2(I-1,1)-C2(I,1))*(C2(I-1,n)-C2(I,2))+C2(I,n);
    end
    if sys.Aligned ==false
        n=3;
        C2=(sys.Columns-C2(I,1))/(C2(I-1,1)-C2(I,1))*(C2(I-1,n)-C2(I,2))+C2(I,n);
    end    
end

end
if sys.Columns>=20
    C2=1;
    if CF.REd>10 && CF.REd<10^2
        C1= 0.8; m= .4;
        
    end
    if CF.REd>10^2 && CF.REd<10^3
        C1=0.683; m=.466;
    end
end

if sys.Aligned == true
    if CF.REd>10 && CF.REd<10^2
        C1= 0.8; m= .4;
        
    end
    if CF.REd>10^2 && CF.REd<10^3
        C1=0.683; m=.466;
    end
    
end
CF.REd
CF.NUD=C1*C2*CF.REd^m*CF.Prf^.36*(CF.Prf/TF.Prf)^(1/4);

CF.hbar=CF.NUD*CF.kf/sys.D;
TF.REd=4*TF.mdot/(sys.N*pi*sys.D*TF.muf);

if TF.REd >= 2300
    TF.xfdh=10*D;
    TF.xfdt=TF.xfdh;
    TF.NUD=0.023*TF.REd^(0.8)*TF.Prf^0.3;
end

if TF.REd<2300
    TF.xfdh=0.05*TF.REd*sys.D;
    TF.xfdt=TF.xfdh*TF.Prf;
    TF.Gzd=sys.D/sys.L*TF.REd*TF.Prf;
    frac1= 3.66/tanh(2.264*TF.Gzd^(-1/3))+1.7*TF.Gzd^(-2/3);
    num1=frac1+0.0499*TF.Gzd*tanh(TF.Gzd^-1);
    den=tanh(2.432*TF.Prf^(1/6)*TF.Gzd^(-1/6));
    TF.NUD=num1/den;
end
TF.hbar=TF.NUD*TF.kf/sys.D;
RCF=1/(CF.hbar*pi*sys.D);
RTF=1/(TF.hbar*pi*sys.D);

Req=RCF+RTF;
TF.L=UA*Req;
CF.L=UA*Req;
end