function res = mtimes(a,b)

if isa(a,'Wavelet_CCIR') == 0
    error('In  A.*B only A can be Wavelet operator');
end

nz = size(b,3);
nt = size(b,4);
    
res = zeros(size(b));

if a.adjoint
    for z = 1:nz
        for k = 1:nt
            res(:,:,z,k) = IWT2_PO(real(b(:,:,z,k)),a.wavScale,a.qmf) + i*IWT2_PO(imag(b(:,:,z,k)),a.wavScale,a.qmf);
        end
    end
else
    for z = 1:nz
        for k = 1:nt
            res(:,:,z,k) = FWT2_PO(real(b(:,:,z,k)),a.wavScale,a.qmf) + i* FWT2_PO(imag(b(:,:,z,k)),a.wavScale,a.qmf);
        end
    end
end


