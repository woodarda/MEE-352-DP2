close all; clear;

Ti=80+273; To=40+273; mdot=[3200/(90*60),3200/(150*60)]; Tbar=(Ti+To)/2;
[Temperature, P,vf, hfg, cpf, muf, kf, Prf]=AW_Interpolation(Tbar);

%-------------------------------------------------------------------------%
%% Conentric Tubes Parameters 
Do=30/1000; % inner tube outer diameter
Di=15/1000; % inner tube inner diameter
D=35/1000;  % outer tube inner diameter
Dair=D+15/1000; % outermost diameter
Dh=D-Do;    % annulus diameter
ks=400;     % convection coefficient
% Surface area of the inner fluid's contact 
% with the tube as a function of L
Ai=pi*Di;
% Surface area of the annulus fluid's contact
% with the tube as a function of L
Ao=pi*Do;   

% number of times that the tube is wraped around 
Nwraps=5;

% changes the fluid that is in the annulus region
% true for design 1 false for design 2
WaterAnnulus=false;

% create a table to track parameters of the geometry
CT=table(Do,Di,D,Dh,ks);
%-------------------------------------------------------------------------%
% NOT USED FOR REPORT

%% Cross Flow Parameters
Rows=4;
Columns=8;
N=Rows*Columns;
L_g=1;
Height=.1;
Width=.1;
ST=Height/Rows;
SL=Width/Columns;
SD=sqrt(SL^2+(ST/2)^2);
A_CS=L_g*Height;
D=ST/2.25;

Aligned=true;
CFP=table(Rows,Columns,N,L_g,Height, Width, ST, SL, SD, A_CS,D,Aligned);
CFP=renamevars(CFP,"A_CS","A");
% NOT USED FOR REPORT
%-------------------------------------------------------------------------%


% create table for the properties of the Chemical
Chem=table(Ti,To,Tbar,mdot,P,vf, hfg, cpf, muf, kf, Prf);
% add the density to the chem table
Chem.Rho=1/Chem.vf;

% calculate the heat transfer rate of the hot fluid
Chem.q=Chem.mdot*Chem.cpf*(Chem.Ti-Chem.To);


Ti=10+273; % water inlet temperature
mdot=0.5; % Water flow rate (kg/s)
Toguess=Chem.To; % initial guess for outlet temp of water
To=[Toguess,Toguess];
% create a table for the parameters of the water
Water=table(Ti, To, mdot);

% use function CalcToH2O to iterate water outlet temperature for 90 mins
% (1)
Water.To(1)=CalcToH2O(Water.To(1),Water.Ti,Water.mdot,Chem.q(1));
% use function CalcToH2O to iterate water outlet temperature for 150 mins
% (2)
Water.To(2)=CalcToH2O(Water.To(2),Water.Ti,Water.mdot,Chem.q(2));
% add tbar_actual to water table 90 min
Water.Tbar(1,1)=(Water.To(1)+Water.Ti)/2;
% add tbar_actual to water table 150 min
Water.Tbar(1,2)=(Water.To(2)+Water.Ti)/2;

% calculate properties of water with actual Tbar for 90 min
[Tbar(1,1), P(1,1),vf(1,1), hfg(1,1), cpf(1,1), muf(1,1), kf(1,1),...
    Prf(1,1)]=AW_Interpolation(Water.Tbar(1,1)); 
% calculate properties of water with actual Tbar for 150 min
[Tbar(1,2), P(1,2),vf(1,2), hfg(1,2), cpf(1,2), muf(1,2), kf(1,2),...
    Prf(1,2)]=AW_Interpolation(Water.Tbar(1,2));

% add properties to the table
Water.P=P; Water.vf=vf;Water.hfg=hfg;Water.cpf=cpf;Water.muf=muf;Water.kf=kf;
Water.Prf=Prf; Water.Rho=1./Water.vf;Water.Tbar=Tbar;

% calculate  heat capacity rate-tios
Cc=Water.mdot*Water.cpf;
Ch=Chem.mdot*Chem.cpf;

% create a table for the 90 min operation
Run90=table(Cc,Ch(1));
Run90=renamevars(Run90,"Var2","Ch"); % rename variable
% create a table for the 150 min operation
Run150=table(Cc,Ch(2)); 
Run150=renamevars(Run150,"Var2","Ch");% rename variable

% add the heat transfer rate into the Run tables
Run90.q=Chem.q;
Run150.q=Chem.q;

% find the Cmins and Cmax using the FindCminCmax function for 90 mins
[Run90.Cmin,Run90.Cmax]=FindCminCmax(Run90.Ch,Run90.Cc(1));
% use Cmin to calculate max heat transfer rate
Run90.qmax=Run90.Cmin*(Chem.Ti-Water.Ti);

% find the Cmins and Cmax using the FindCminCmax function for 150 mins
[Run150.Cmin,Run150.Cmax]=FindCminCmax(Run150.Ch,Run150.Cc(1));
% use Cmin to calculate max heat transfer rate
Run150.qmax=Run150.Cmin*(Chem.Ti-Water.Ti);

% calculate 90 min efficiency
Eff(1,1)=Run90.q(1)/Run90.qmax;
% calculate 150 min efficiency
Eff(1,2)=Run150.q(2)/Run150.qmax;

% create variable in table for efficiency 
Run90.Eff=Eff(1);
Run150.Eff=Eff(2);

% Calculate the heate transfer rate-tios
Run90.Cr=Run90.Cmin/Run90.Cmax;
Run150.Cr=Run150.Cmin/Run150.Cmax;

% Calculate the NTU for the parallel flows and and counter flow of 90 min
Run90.NTUpar=-log(1-Run90.Eff*(1+Run90.Cr))/(1+Run90.Cr);
Run90.NTUcntr=log((Run90.Eff-1)/(Run90.Eff*Run90.Cr-1))/(Run90.Cr-1);

% Calculate the NTU for the parallel flows and and counter flow of 150 min
Run150.NTUpar=-log(1-Run150.Eff*(1+Run150.Cr))/(1+Run150.Cr);
Run150.NTUcntr=log((Run150.Eff-1)/(Run150.Eff*Run150.Cr-1))/(Run90.Cr-1);

% NTUpar is complex so desired heat transfer is not possible 
% with that exchanger.

% Calculate UA for the counter flow of 90 min run
Run90.UAcntr=Run90.NTUcntr*Run90.Cmin;

% Calculate UA for the parallel and counter flow of 150 min run
Run150.UAcntr=Run150.NTUcntr*Run150.Cmin;
Run150.UApar=Run150.NTUpar*Run90.Cmin;

% NOT USED FOR REPORT DESIGNS
%-------------------------------------------------------------------------%
Run90.NTUCFcmaxMix= -log(1+log(1-Run90.Eff*Run90.Cr)/Run90.Cr);
Run90.NTUCFcmaxUnmix= -log(Run90.Cr*log(1-Run90.Eff)+1)/Run90.Cr;

Run150.NTUCFcmaxMix= -log(1+log(1-Run150.Eff*Run150.Cr)/Run150.Cr);
Run150.NTUCFcmaxUnmix= -log(Run150.Cr*log(1-Run150.Eff)+1)/Run150.Cr;


Run90.UACFcmaxUnmix=Run90.NTUCFcmaxUnmix*Run90.Cmin;

Run150.UACFcmaxMix=Run150.NTUCFcmaxMix*Run150.Cmin;
Run150.UACFcmaxUnmix=Run150.NTUCFcmaxUnmix*Run90.Cmin;
% NOT USED FOR REPORT DESIGNS
%-------------------------------------------------------------------------%

% UA for Concentric tubes with counter flow was smallest for 90 min run
% Calculate L using the Concentric tube function using the 90 min water
% properties and chemical mass flow rate
if WaterAnnulus==true
    [L,Annulus,Circle]=ConcentricTubes(Run150.UAcntr,Water.muf(1),...
        Water.kf(1),Water.Prf(1),Water.To(1),Water.Ti,Chem.muf,Chem.kf,...
        Chem.Prf,Chem.To,Chem.Ti,CT,Water.mdot,Chem.mdot(1),ks);
HAd="------------------------"
%     [L,Annulus,Circle]=ConcentricTubes(Run150.UAcntr,Water.muf(2),...
%         Water.kf(2),Water.Prf(2),Water.To(2),Water.Ti,Chem.muf,Chem.kf,...
%         Chem.Prf,Chem.To,Chem.Ti,CT,Water.mdot,Chem.mdot(2),ks);
end
% Calculate L using the Concentric tube function using the 150 min water
% properties and chemical mass flow rate
if WaterAnnulus==false
    [L,Annulus,Circle]=ConcentricTubes(Run90.UAcntr,Chem.muf,Chem.kf,...
        Chem.Prf,Chem.To,Chem.Ti,Water.muf(2),Water.kf(2),Water.Prf(2),...
        Water.To(2),Water.Ti,CT,Chem.mdot(2),Water.mdot,ks);
end
% length of the pipe needed is stored in annulus and circle table
Annulus.L=L;
Circle.L=L;

% NOT USED FOR REPORT DESIGNS
%-------------------------------------------------------------------------%
delta=100;
Z=0;
Tol=5;
while Z==0
if CFP.Aligned==true
    [CF,TF]=CrossFlow(Run150.UACFcmaxMix,Water,Chem,CFP); 
end

Z=TF.L<CFP.L_g*(1+Tol) && TF.L>CFP.L_g*(1-Tol);
CFP.L_g=TF.L;
CFP.A=L_g*Height;
end
% NOT USED FOR REPORT DESIGNS
%-------------------------------------------------------------------------%

% the distance that each pipe's center will be as it wraps around
dx=Dair/1.8; 
% The length of the heat exchanger
LCT=Annulus.L/Nwraps;
% width has Nwraps-1 gaps + the outer most diameter times the number of
% wraps
WidthCT=(Nwraps-1)*(dx-Dair)+Dair*Nwraps;

% NOT USED FOR REPORT DESIGNS
%-------------------------------------------------------------------------%
% Create a table to view the important parameters for comparison
varname={'Operation Time (mins)','Operation Cost (per day)', 'Length (m)',...
    'Cross Section Area (m^2)','Floor Footprint (m^2)','Volume (m^3)',...
    'Chemical Temperature Inlet (^o C)','Chemical Temperature Outlet (^o C)',...
    'Water Temperature Inlet (^o C)','Water Temperatrue Outlet (^o C)'};
rowname={'Concentric Tubes','Cross Flow Bank of Tubes'};
var1=string([90;150]);
var2=string(['$1,000';'$2,000']);
var3=[LCT;CF.L];
var4=[WidthCT*Dair;Height*Width];
var5=[LCT*WidthCT;CF.L*Width];
var6=[LCT*Dair*WidthCT;CF.L*Width*Height];
var7=[Chem.Ti;Chem.Ti];
var8=[Chem.To;Chem.To];
var9=[Water.Ti;Water.Ti];
var10=[Water.To(1);Water.To(2)];
% NOT USED FOR REPORT DESIGNS
%-------------------------------------------------------------------------%


