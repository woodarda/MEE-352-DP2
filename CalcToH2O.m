function [To_des] = CalcToH2O(To_des,Ti_des,mdot_des,q)
%CALCTOH2O Calculates the outlet temperature and iterates
%   This takes in the properties of the chem and the guessed water temp and
%   calculates the water outlet temperature using qh=qc. It will then
%   iterate until it is within the desired range of Tol and returns the
%   outlet temperature of the water

% when z=0 it will iterate. the calculated tbar_actual (Tbara) is within
% the desired range it changes to true
Z=0;
Tol=.05/100; % desired tolerance for the Tbara
while Z==0
Tbar=(To_des+Ti_des)*0.5; % calculate guess for average temperature
[Temp, P,vf, hfg,cpf_des,muf,kf,Prf]=AW_Interpolation(Tbar); %#ok<ASGLU> 

To_des=q/(mdot_des*cpf_des)+Ti_des; % solve qh = qc for To water

Tbara=To_des/2+Ti_des/2; % calculate Tbar_act

Z=Tbara<Tbar*(1+Tol) && Tbara>Tbar*(1-Tol); % test if within tolerance 

end

end

