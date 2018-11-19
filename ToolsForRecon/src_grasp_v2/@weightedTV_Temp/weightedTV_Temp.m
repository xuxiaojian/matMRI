function  res = weightedTV_Temp(weights)

%res = TV_Temp()
%
% Implements a difference operator along the time dimensin for dynamic MRI
%
% Ricardo Otazo 2008
%
% Cihat Eldeniz modified TV to weighted TV [May 2017]

res.adjoint = 0;
res.w = weights; % wTV
res = class(res,'weightedTV_Temp'); % wTV

