function actualToAnalyze(data,filename,x,y,z)

[width,height,numberOfFrames,numberOfSlices] = size(data);

data = ipermute(data,[1 2 4 3]);
data = flip(data,1);
data = reshape(data,[width height numberOfSlices*numberOfFrames]);
data = ipermute(data,[2 3 1]);

analyzeWriter(data,x,y,z,16,filename);