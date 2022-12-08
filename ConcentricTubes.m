function [L,Annulus,Circle] = ConcentricTubes(UA,Amuf,Akf,...
    APrf,ATo,ATi,Cmuf,Ckf,CPrf,CTo,CTi,Tube,mdotA,mdotC,ks)
%CONCENTRICTUBES Summary of this function goes here
%   Detailed explanation goes here
Annulus=table(Amuf,Akf,APrf,ATo,ATi);
Annulus=renamevars(Annulus,["Amuf","Akf","APrf","ATo","ATi"],["muf","kf",...
    "Prf","To","Ti"]);

Circle=table(Cmuf,Ckf,CPrf,CTo,CTi);
Circle=renamevars(Circle,["Cmuf","Ckf","CPrf","CTo","CTi"],["muf","kf",...
    "Prf","To","Ti"]);

Annulus.REd=4*mdotA/(Tube.Dh*pi*Annulus.muf);
Annulus.xfdh=10*Tube.Dh;
Annulus.xfdt=Annulus.xfdh;

Circle.REd=4*mdotC/(Tube.Di*pi*Circle.muf);
Circle.xfdh=10*Tube.Di;
Circle.xfdt=Circle.xfdh;
% Annulus.f=(0.79*log(Annulus.REd)-1.64)^-2;
% Circle.f=(0.79*log(Annulus.REd)-1.64)^-2;
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

% Annulus.NUdbar=(Annulus.f/8*(Annulus.REd-1000)*Annulus.Prf)/...
%     (1+12.7*(Annulus.f/8)^(1.2)*Annulus.Prf^(2/3)-1);
% Circle.NUdbar=(Circle.f/8*(Circle.REd-1000)*Circle.Prf)/...
%     (1+12.7*(Circle.f/8)^(1.2)*Circle.Prf^(2/3)-1);


Annulus.hbar=Annulus.NUdbar*Annulus.kf/Tube.Di;
Circle.hbar=Circle.NUdbar*Circle.kf/Circle.muf;
Rw=log(Tube.Do/Tube.Di)/(2*pi*ks);
RhA=1/(Annulus.hbar*pi*Tube.Do);
RhC=1/(Circle.hbar*pi*Tube.Di);
Annulus.Req=RhA+RhC+Rw;

L=UA*Annulus.Req;
Annulus.UA=UA;
Circle.UA=UA;
Circle.Req=Annulus.Req;
end

