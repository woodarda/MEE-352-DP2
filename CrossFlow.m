function [CF,TF] = CrossFlow(UA,CF,TF,sys)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
CF.vf=CF.vf(2); CF.P=CF.P(2); CF.hfg=CF.hfg(2); CF.cpf=CF.cpf(2);
CF.muf=CF.muf(2);CF.kf=CF.kf(2);CF.Prf=CF.Prf(2);CF.Rho=CF.Rho(2);

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
end

if sys.Aligned == true
    if CF.REd>10 && CF.REd<10^2
        C1= 0.8; m= .4;
        
    end
    if CF.REd>10^2 && CF.REd<10^3
        C1=0.683; m=.466;
    end
    
end

CF.NUD=C1*C2*CF.REd^m*CF.Prf^.36*(CF.Prf/TF.Prf)^(1/4);

CF.hbar=CF.NUD*CF.kf/sys.D;

TF = 0;
end