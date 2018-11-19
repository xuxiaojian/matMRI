function res = mtimes(a,u)

[~,Nx,Np] = size(u);
    
if a.adjoint
    res = cat(1,zeros(1,Nx,Np),u(1:end-1,:,:)) - cat(1,u(1:end-1,:,:),zeros(1,Nx,Np));
else
    res = cat(1,u(2:end,:,:),u(end,:,:)) - u;
end