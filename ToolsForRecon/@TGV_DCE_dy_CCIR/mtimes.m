function res = mtimes(a,u)

[~,Nx,Np,Nc] = size(u);
    
if a.adjoint
    res = cat(1,zeros(1,Nx,Np,Nc),u(1:end-1,:,:,:)) - cat(1,u(1:end-1,:,:,:),zeros(1,Nx,Np,Nc));
else
    res = cat(1,u(2:end,:,:,:),u(end,:,:,:)) - u;
end