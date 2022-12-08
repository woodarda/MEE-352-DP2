function [To_des] = CalcToH2O(To_des,Ti_des,mdot_des,q)
%CALCTOH2O Summary of this function goes here
%   Detailed explanation goes here
Z=0;
Tol=.05/100;
while Z==0
Tbar=(To_des+Ti_des)*0.5;
[Temp, P,vf, hfg,cpf_des,muf,kf,Prf]=AW_Interpolation(Tbar); %#ok<ASGLU> 

To_des=q/(mdot_des*cpf_des)+Ti_des;

Tbara=To_des/2+Ti_des/2;

Z=Tbara<Tbar*(1+Tol) && Tbara>Tbar*(1-Tol);

end

end

