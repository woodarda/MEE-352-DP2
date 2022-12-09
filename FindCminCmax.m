function [Cmin,Cmax] = FindCminCmax(Ch,Cc)
%FindCminCmax calculates the heat capacity rate-tios 
%   If the ot fluid's heat capacity rate-tio is lower it is Cmin and the 
% cold fluid is Cmax. if the cold fluid has a lower C it is the Cmin and
% the hot fluid is Cmax
if Ch<Cc
    Cmin=Ch;
    Cmax=Cc;
end

if Cc<Ch
    Cmin=Cc;
    Cmax=Ch;
end

end

