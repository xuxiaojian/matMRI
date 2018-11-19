function ress = mtimes(a,bb)

if a.adjoint,   
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point    
    res = cell(length(bb),1);    
    ress = zeros([a.imSize length(bb)]);  
    for tt=1:length(bb),
        for ch=1:size(a.b1,3),
            b = bb{tt}(:,:,ch).*sqrt(a.w{tt});
            res{tt}(:,:,ch) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize)),a.imSize(1),a.imSize(2));                                    
        end    
        % compensate for undersampling factor
        res{tt}=res{tt}*size(a.b1,1)*pi/2/size(a.w{tt},2);
        
        % coil combination for each time point              
        for y = 1:size(res{tt},1)
            for x = 1:size(res{tt},2)
                if abs(a.b1(y,x))>0
                    ress(y,x,tt)=sum(res{tt}(y,x,:).*conj(a.b1(y,x,:)),3)./sum(abs((a.b1(y,x,:))).^2,3); 
                end
            end
        end
    end         
else
    % Cartesian image to multicoil non-Cartesian k-space
    ress = cell(size(bb,3),1);
    for tt=1:size(bb,3),
        ress{tt} = zeros([a.dataSize(1) a.dataSize(2) size(a.b1,3)]);
        for ch=1:size(a.b1,3),
            res=bb(:,:,tt).*a.b1(:,:,ch);
            ress{tt}(:,:,ch) = reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize)),a.dataSize(1),a.dataSize(2)).*sqrt(a.w{tt});
        end
    end
end

