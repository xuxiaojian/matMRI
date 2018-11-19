function K = GenerateRadialSpokesNew(nSpokes, ndata, spokes_ind)

% M_PI = 3.14159265358979323846;

K = zeros(ndata, nSpokes);

R = 0.5;

baseRes = 320; % mm.
dk = 1/baseRes;

k = -R:dk/2:R-dk/2; % Oversampled
% Set the DC point to absolute zero [not to eps or anything else]
k(ndata/2+1) = 0;

for i = 1:nSpokes

    A = inc_golden(spokes_ind(i) - 1);
     
%     kx = R*[cos(A), cos(A + M_PI)];
%     ky = R*[sin(A), sin(A + M_PI)];
%     
%     kx1 = kx(1);
%     kx2 = kx(2);
% 
%     xi = linspace(kx1, kx2, ndata);
%     yi = linspace(ky(1), ky(2), ndata);
%    
%     K(:, i) = xi(:) + 1i * yi(:);

    kx = k*cos(A);
    ky = k*sin(A);
    
    K(:, i) = kx(:) + 1i * ky(:);
end