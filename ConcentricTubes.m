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

Annulus.hbar=Annulus.NUdbar*Annulus.kf/Annulus.muf;
Circle.hbar=Circle.NUdbar*Circle.kf/Circle.muf;
Rw=log(Tube.Do/Tube.Di)/(2*pi*ks);
RhA=1/(Annulus.hbar*pi*Tube.Do);
RhC=1/(Circle.hbar*pi*Tube.Di);
Req=RhA+RhC+Rw;

L=UA*Req;

end

