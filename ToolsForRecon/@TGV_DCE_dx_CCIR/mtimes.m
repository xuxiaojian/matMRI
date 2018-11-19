function res = mtimes(a,u)

[Ny,~,Np,Nc] = size(u);
    
if a.adjoint
    res = cat(2,zeros(Ny,1,Np,Nc),u(:,1:end-1,:,:)) - cat(2,u(:,1:end-1,:,:),zeros(Ny,1,Np,Nc));
else
    res = cat(2,u(:,2:end,:,:),u(:,end,:,:)) - u;
end


