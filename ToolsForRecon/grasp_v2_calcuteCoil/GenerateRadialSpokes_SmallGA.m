function K = GenerateRadialSpokes_SmallGA(nSpokes, ndata)

M_PI = 3.14159265358979323846;

K = zeros(ndata, nSpokes);

R = 0.5;

for i = 1:nSpokes

%     A = inc_golden(spokes_ind(i) - 1);    
%     A = pi*i/nSpokes*cos(pi*i); % cosine
%     A = pi*i/nSpokes; % ramp
    A = pi*68.753882025/180*(i-1); % small GA
        
    kx = R*[cos(A), cos(A + M_PI)];
    ky = R*[sin(A), sin(A + M_PI)];
    
    kx1 = kx(1);
    kx2 = kx(2);

    xi = linspace(kx1, kx2, ndata);
    yi = linspace(ky(1), ky(2), ndata);
   
    K(:, i) = xi(:) + 1i * yi(:);
    
end