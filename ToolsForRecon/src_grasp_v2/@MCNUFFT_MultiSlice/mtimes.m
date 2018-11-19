function ress = mtimes(a,bb)

% Multislice: The sensitivity matrix a.b1 is Ny * Nx * Nc * Nz


if a.adjoint,   
    % Multislice: The k-space data bb is Ndata * NspokesPerPhase * Nc * Nz * Nt
    
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    res = zeros([size(a.b1) size(bb,5)]); % Multislice
    for z = 1:size(bb,4) % Multislice
        for tt=1:size(bb,5),
            for ch=1:size(bb,3),
                b = bb(:,:,ch,z,tt).*a.w(:,:,tt);
                res(:,:,ch,z,tt) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize)),a.imSize(1),a.imSize(2));                                    
            end        
        end 
    end
    % compensate for undersampling factor
    res=res*size(a.b1,1)*pi/2/size(a.w,2);
    % coil combination for each time point
%     ress = zeros([size(res,1) size(res,2) size(bb,4)]);
%     for tt=1:size(bb,4),
%         for y = 1:size(res,1)
%             for x = 1:size(res,2)
%                 if abs(a.b1(y,x))>0
%                     ress(y,x,tt)=sum(res(y,x,:,tt).*conj(a.b1(y,x,:)),3)./sum(abs((a.b1(y,x,:))).^2,3); 
%                 end
%             end
%         end
%     end    
    ress = squeeze(sum(res.*repmat(conj(a.b1),[1 1 1 1 size(bb,5)]),3))./repmat(squeeze(sum(abs(a.b1).^2,3)),[1 1 1 size(bb,5)]); % Multislice ==> ress = Ny * Nx * Nz * Nt
else
    % Multislice: The image-space data bb is Ny * Nx * Nz * Nt
    
    % Cartesian image to multicoil non-Cartesian k-space
    ress = zeros([a.dataSize(1) a.dataSize(2) size(a.b1,3) size(bb,3) size(bb,4)]);
    for z = 1:size(bb,3) % Multislice
        for tt=1:size(bb,4),
            for ch=1:size(a.b1,3),
                res=bb(:,:,z,tt).*a.b1(:,:,ch,z);
                ress(:,:,ch,z,tt) = reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);
            end
        end
    end
end

