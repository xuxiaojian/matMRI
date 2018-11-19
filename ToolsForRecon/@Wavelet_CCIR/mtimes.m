function res = mtimes(a,b)

if isa(a,'Wavelet_CCIR') == 0
    error('In  A.*B only A can be Wavelet operator');
end

nt = size(b,3);
    
res = zeros(size(b));

if a.adjoint
    for k = 1:nt
        res(:,:,k) = IWT2_PO(real(b(:,:,k)),a.wavScale,a.qmf) + i*IWT2_PO(imag(b(:,:,k)),a.wavScale,a.qmf);
    end
else
    for k = 1:nt
        res(:,:,k) = FWT2_PO(real(b(:,:,k)),a.wavScale,a.qmf) + i* FWT2_PO(imag(b(:,:,k)),a.wavScale,a.qmf);
    end
end


