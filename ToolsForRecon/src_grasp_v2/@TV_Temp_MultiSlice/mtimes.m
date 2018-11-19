function res = mtimes(a,b)

% Cihat: b is 4D with dimensions Ny * Nx * Nz * Nt

if a.adjoint
    res = adjDz(b);
else
    res = b(:,:,:,[2:end,end]) - b;
end

function y = adjDz(x)
y= x(:,:,:,[1,1:end-1]) - x;
y(:,:,:,1) = -x(:,:,:,1);
y(:,:,:,end) = x(:,:,:,end-1);