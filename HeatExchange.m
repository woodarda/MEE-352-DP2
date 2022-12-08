close all; clear; clc;

Ti=80+273; To=40+273; mdot=[3200/(90*60),3200/(150*60)]; Tbar=(Ti+To)/2;
[Temperature, P,vf, hfg, cpf, muf, kf, Prf]=AW_Interpolation(Tbar);

%-------------------------------------------------------------------------%
%% Conentric Tubes Parameters 
Do=50/1000;
Di=45/1000;
D=65/1000;
Dh=D-Do;
ks=400;
Ai=pi*Di;
Ao=pi*Do;


WaterAnnulus=true;

CT=table(Do,Di,D,Dh,ks);
%-------------------------------------------------------------------------%
%% Cross Flow Parameters
Rows=6;
Columns=14;
N=Rows*Columns;
L_g=5;
Height=10;
Width=6;
ST=Height/Rows;
SL=Width/Columns;
SD=sqrt(SL^2+(ST/2)^2);
A_CS=L_g*Height;
D=ST/2.5;

Aligned=true;
CFP=table(Rows,Columns,N,L_g,Height, Width, ST, SL, SD, A_CS,D,Aligned);
CFP=renamevars(CFP,"A_CS","A");
%-------------------------------------------------------------------------%
Chem=table(Ti,To,Tbar,mdot,P,vf, hfg, cpf, muf, kf, Prf);
Chem.Rho=1/Chem.vf;

Chem.q=Chem.mdot*Chem.cpf*(Chem.Ti-Chem.To);

Ti=10+273; % water inlet temperature
mdot=0.5; % Water flow rate (kg/s)
Toguess=Chem.To; % initial guess for outlet temp of water
To=[Toguess,Toguess];
Water=table(Ti, To, mdot);

Water.To(1)=CalcToH2O(Water.To(1),Water.Ti,Water.mdot,Chem.q(1));
Water.To(2)=CalcToH2O(Water.To(2),Water.Ti,Water.mdot,Chem.q(2));
Water.Tbar(1,1)=(Water.To(1)+Water.Ti)/2;
Water.Tbar(1,2)=(Water.To(2)+Water.Ti)/2;

[Tbar(1,1), P(1,1),vf(1,1), hfg(1,1), cpf(1,1), muf(1,1), kf(1,1),...
    Prf(1,1)]=AW_Interpolation(Water.Tbar(1,1));  
[Tbar(1,2), P(1,2),vf(1,2), hfg(1,2), cpf(1,2), muf(1,2), kf(1,2),...
    Prf(1,2)]=AW_Interpolation(Water.Tbar(1,2));

Water.P=P; Water.vf=vf;Water.hfg=hfg;Water.cpf=cpf;Water.muf=muf;Water.kf=kf;
Water.Prf=Prf; Water.Rho=1./Water.vf;Water.Tbar=Tbar;
Cc=Water.mdot*Water.cpf;
Ch=Chem.mdot*Chem.cpf;

Run90=table(Cc,Ch(1));
Run90=renamevars(Run90,"Var2","Ch");
Run150=table(Cc,Ch(2));
Run150=renamevars(Run150,"Var2","Ch");

Run90.q=Chem.q;
Run150.q=Chem.q;

[Run90.Cmin,Run90.Cmax]=FindCminCmax(Run90.Ch,Run90.Cc(1));
Run90.qmax=Run90.Cmin*(Chem.Ti-Water.Ti);

[Run150.Cmin,Run150.Cmax]=FindCminCmax(Run150.Ch,Run150.Cc(1));
Run150.qmax=Run150.Cmin*(Chem.Ti-Water.Ti);

Eff(1,1)=Run90.q(1)/Run90.qmax;
Eff(1,2)=Run150.q(2)/Run150.qmax;

Run90.Eff=Eff(1);
Run150.Eff=Eff(2);

Run90.Cr=Run90.Cmin/Run90.Cmax;
Run150.Cr=Run150.Cmin/Run150.Cmax;

Run90.NTUpar=-log(1-Run90.Eff*(1+Run90.Cr))/(1+Run90.Cr);
Run90.NTUcntr=log((Run90.Eff-1)/(Run90.Eff*Run90.Cr-1))/(Run90.Cr-1);

Run150.NTUpar=-log(1-Run150.Eff*(1+Run150.Cr))/(1+Run150.Cr);
Run150.NTUcntr=log((Run150.Eff-1)/(Run150.Eff*Run150.Cr-1))/(Run90.Cr-1);

Run90.UAcntr=Run90.NTUcntr*Run90.Cmin;

Run150.UAcntr=Run150.NTUcntr*Run150.Cmin;
Run150.UApar=Run150.NTUpar*Run90.Cmin;

Run90.NTUCFcmaxMix= -log(1+log(1-Run90.Eff*Run90.Cr)/Run90.Cr);
Run90.NTUCFcmaxUnmix= -log(Run90.Cr*log(1-Run90.Eff)+1)/Run90.Cr;

Run150.NTUCFcmaxMix= -log(1+log(1-Run150.Eff*Run150.Cr)/Run150.Cr);
Run150.NTUCFcmaxUnmix= -log(Run150.Cr*log(1-Run150.Eff)+1)/Run150.Cr;


Run90.UACFcmaxUnmix=Run90.NTUCFcmaxUnmix*Run90.Cmin;

Run150.UACFcmaxMix=Run150.NTUCFcmaxMix*Run150.Cmin;
Run150.UACFcmaxUnmix=Run150.NTUCFcmaxUnmix*Run90.Cmin;

% UA for Concentric tubes with counter flow was smallest for 90 min run
if WaterAnnulus==true
    [L,Annulus,Circle]=ConcentricTubes(Run90.UAcntr,Water.muf(1),...
        Water.kf(1),Water.Prf(1),Water.To(1),Water.Ti,Chem.muf,Chem.kf,...
        Chem.Prf,Chem.To,Chem.Ti,CT,Water.mdot,Chem.mdot(1),ks);
end

if WaterAnnulus==false
    [L,Annulus,Circle]=ConcentricTubes(Run90.UAcntr,Chem.muf,Chem.kf,...
        Chem.Prf,Chem.To,Chem.Ti,Water.muf(1),Water.kf(1),Water.Prf(1),...
        Water.To(1),Water.Ti,CT,Chem.mdot(1),Water.mdot,ks);
end
Annulus.L=L;
Circle.L=L;
delta=100;

while delta>5

    if CFP.Aligned==true
        [CF,TF]=CrossFlow(Run150.UACFcmaxMix,Water,Chem,CFP); 
    end


L_CF=TF.L/N;
TF_RE=TF.REd;
CF_RE=TF.REd;
delta=abs(L_CF-CFP.L_g)/CFP.L_g;
CFP.L_g=L_CF;
CFP.A=CFP.L_g*Height;
end

TF.LN=TF.L/N;
