function [L,Annulus,Circle] = ConcentricTubes(UA,Amuf,Akf,...
    APrf,ATo,ATi,Cmuf,Ckf,CPrf,CTo,CTi,Tube,mdotA,mdotC,ks)
%CONCENTRICTUBES Calculates the L for the Concentric tubes heat exchanger
% Calculates the heat transfer rates for the concentric tube circular
% region and the annulus region as well as the conduction restistance
% between them. Flow is Turbulent and entrance lengths are assumed to be
% negligible. Annulus region can be treated as circular with Dh=D-Do

% create a table for the properties of the annulus region fluid
Annulus=table(Amuf,Akf,APrf,ATo,ATi);
Annulus=renamevars(Annulus,["Amuf","Akf","APrf","ATo","ATi"],["muf","kf",...
    "Prf","To","Ti"]);

% creates a table for the properties of the circular fluid flow
Circle=table(Cmuf,Ckf,CPrf,CTo,CTi);
Circle=renamevars(Circle,["Cmuf","Ckf","CPrf","CTo","CTi"],["muf","kf",...
    "Prf","To","Ti"]);

% Calculate the Reynolds number and the entry length effects to later
% verify they are small.
Annulus.REd=4*mdotA/(Tube.Dh*pi*Annulus.muf);
Annulus.xfdh=10*Tube.Dh;
Annulus.xfdt=Annulus.xfdh;

% Calculate the Reynolds number and the entry length effects to later
% verify they are small.
Circle.REd=4*mdotC/(Tube.Di*pi*Circle.muf);
Circle.xfdh=10*Tube.Di;
Circle.xfdt=Circle.xfdh;

% Checks the temperature for which is hot and cold so that the correct n
% can be used to calculate the Nusselt number. uses Turbulent calculation
% for NUd
if Annulus.Ti<Circle.Ti
    n=.4; 
    Annulus.NUdbar=0.023*Annulus.REd^(0.8)*Annulus.Prf^(n);
    n=.3;
    Circle.NUdbar=0.023*Circle.REd^(0.8)*Circle.Prf^(n);
end
if Annulus.Ti>Circle.Ti
    n=.3;
    Annulus.NUdbar=0.023*Annulus.REd^(0.8)*Annulus.Prf^(n);
    n=0.4;
    Circle.NUdbar=0.023*Circle.REd^(0.8)*Circle.Prf^(n);
end

% Calculate the convection coefficents for the Annulus and Circular region
Annulus.hbar=Annulus.NUdbar*Annulus.kf/Tube.Di;
Circle.hbar=Circle.NUdbar*Circle.kf/Circle.muf;
% calculate the conduction resistance
Rw=log(Tube.Do/Tube.Di)/(2*pi*ks)
% Calculate the annulus convection resistance
RhA=1/(Annulus.hbar*pi*Tube.Do);
% Calculate the circular region convection resistance
RhC=1/(Circle.hbar*pi*Tube.Di);
% Calculate the equivalent resistance 
Annulus.Req=RhA+RhC+Rw;

% Calculate L using UA=1/L*Req'
L=UA*Annulus.Req;
Annulus.UA=UA;
Circle.UA=UA;
Circle.Req=Annulus.Req;
end

