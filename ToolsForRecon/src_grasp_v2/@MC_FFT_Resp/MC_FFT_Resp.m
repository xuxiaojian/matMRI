function  res = MC_FFT_Resp(b1,missingIndices)
% b1 keeps single-slice coil sensitivity maps for 32 channels [Ny*Nx*Nc]
% s is the index of the slice of interest.
res.adjoint = 0;
res.b1 = b1;
res.missingIndices = missingIndices;
res = class(res,'MC_FFT_Resp');

