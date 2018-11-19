function ress = mtimes_withScaling(a,bb)

if a.adjoint,    
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    res = zeros([size(a.b1) size(bb,4)]);
    
    if ~evalin('base','exist(''scalingFactors'',''var'')')
        scalingFactors = zeros(size(bb,3),1);
        bScalingFactorsExist = false;
    else        
        bScalingFactorsExist = true;
    end    
    
    for tt=1:size(bb,4),
        for ch=1:size(bb,3),
            b = bb(:,:,ch,tt).*a.w(:,:,tt);
            res(:,:,ch,tt) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize)),a.imSize(1),a.imSize(2));            
            if tt==1 && ~bScalingFactorsExist
                scalingFactors(ch) = max(max(abs(res(:,:,ch,tt))));
            end            
        end        
    end  
    if ~bScalingFactorsExist        
        assignin('base','scalingFactors',scalingFactors);
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
    ress = squeeze(sum(res.*repmat(conj(a.b1),[1 1 1 size(bb,4)]),3))./repmat(sum(abs(a.b1).^2,3),[1 1 size(bb,4)]);
else
    % Cartesian image to multicoil non-Cartesian k-space
    ress = zeros([a.dataSize(1) a.dataSize(2) size(a.b1,3) size(bb,3)]);
    for tt=1:size(bb,3),
        for ch=1:size(a.b1,3),
            res=bb(:,:,tt).*a.b1(:,:,ch);
            ress(:,:,ch,tt) = reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);
        end
    end
end

