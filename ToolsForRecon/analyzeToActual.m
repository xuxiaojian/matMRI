function [data,x,y,z] = analyzeToActual(filename)

hdr = analyze75info(filename);
data = double(analyze75read(filename));
x = double(hdr.PixelDimensions(1));
y = double(hdr.PixelDimensions(2));
z = double(hdr.PixelDimensions(3));

[height,width,numberOfSlices] = size(data);
data = reshape(data,[height width numberOfSlices 1]);
data = flipdim(data,1);
data = permute(data,[1 2 4 3]); 