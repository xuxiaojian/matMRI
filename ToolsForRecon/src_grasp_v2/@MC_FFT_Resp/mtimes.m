function ress = mtimes(a,bb)

if a.adjoint
    % bb is k-space data of size Nx*Ny*Nz*Nc*Np where Nc is the number of
    % coils and Np is the number of phases.
    [Nx,Ny,Nz,Nc,Np] = size(bb);
    ress = zeros(Ny,Nx,Nz,Np);
    for ph=1:Np
        for ch=1:Nc
            I = abs(fftshift(ifft(ifftshift(fftshift(ifft(ifftshift(fftshift(ifft(ifftshift(bb(:,:,:,ch,ph),1),[],1),1),2),[],2),2),3),[],3),3));
            I = permute(flip(flip(I,1),3),[2 1 3]);
            ress(:,:,:,ph) = sqrt(ress(:,:,:,ph).^2+I.^2);
        end
    end
else
    % bb is multi-phase image-space data of size Ny*Nx*Nz*Np where Np is
    % the number of phases.
    [Ny,Nx,Nz,Np] = size(bb);    
    % The kx*ky plane is transposed with respect to the image space. So, we
    % use Nx,Ny,... rather than Ny,Nx,...
    ress = zeros(Nx,Ny,Nz,size(a.b1,4),Np);
    for ph=1:Np
        for ch=1:size(a.b1,4)
            res = bb(:,:,:,ph).*a.b1(:,:,:,ch);
            res = ipermute(flip(flip(res,1),3),[2 1 3]);
            res = fftshift(ifft(ifftshift(res,1),[],1),1);
            res = fftshift(ifft(ifftshift(res,2),[],2),2);
            res = fftshift(ifft(ifftshift(res,3),[],3),3);
            res(:,:,a.missingIndices(a.missingIndices(:,1)==ph,2)) = 0;
            ress(:,:,:,ch,ph) = res;
            clear res
        end        
    end
end

