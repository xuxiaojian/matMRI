function ress = mtimes(a,bb)

if a.adjoint,    
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    res = zeros([size(a.b1) size(bb,4) size(bb,5)]);    
        
    for c = 1:size(bb,5)
        for tt=1:size(bb,4),
            for cha=1:size(bb,3),
                b = bb(:,:,cha,tt,c).*a.w(:,:,tt,c);
                res(:,:,cha,tt,c) = reshape(nufft_adj(b(:),a.st{tt,c})/sqrt(prod(a.imSize)),a.imSize(1),a.imSize(2));                                    
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
    ress = squeeze(sum(res.*repmat(conj(a.b1),[1 1 1 size(bb,4) size(bb,5)]),3))./repmat(sum(abs(a.b1).^2,3),[1 1 size(bb,4) size(bb,5)]);
else
    % Cartesian image to multicoil non-Cartesian k-space
    ress = zeros([a.dataSize(1) a.dataSize(2) size(a.b1,3) size(bb,3) size(bb,4)]);
    for c = 1:size(bb,4)
        for tt=1:size(bb,3),
            for cha=1:size(a.b1,3),
                res=bb(:,:,tt,c).*a.b1(:,:,cha);
                ress(:,:,cha,tt,c) = reshape(nufft(res,a.st{tt,c})/sqrt(prod(a.imSize)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt,c);
            end
        end
    end
end

