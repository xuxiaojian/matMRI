%  Save Analyze7.5 header and image file

function analyzeWriter(data, pixdim_x, pixdim_y, pixdim_z, data_type, filename)

%%%%%%%%%%%%%%%%%%%fabricate the header

hdr.header_key.sizeof_hdr = 348;
hdr.header_key.data_type = zeros(10, 1, 'int8');
hdr.header_key.db_name = zeros(18, 1, 'int8');
hdr.header_key.extents = 0; 
hdr.header_key.session_error = 0;
hdr.header_key.regular = 0;
hdr.header_key.hkey_un0 = 0;

hdr.image_dimension.dim(1) = 4;
hdr.image_dimension.dim(2) = size(data, 3);
hdr.image_dimension.dim(3) = size(data, 2);
hdr.image_dimension.dim(4) = size(data, 1);
hdr.image_dimension.dim(5) = 1;
hdr.image_dimension.dim(6) = 0;
hdr.image_dimension.dim(7) = 0;
hdr.image_dimension.dim(8) = 0;

hdr.image_dimension.unused8 = 0;
hdr.image_dimension.unused9 = 0;
hdr.image_dimension.unused10 = 0;
hdr.image_dimension.unused11 = 0;
hdr.image_dimension.unused12 = 0;
hdr.image_dimension.unused13 = 0;
hdr.image_dimension.unused14 = 0;

hdr.image_dimension.datatype = data_type;

switch hdr.image_dimension.datatype
   case   1,
      hdr.image_dimension.bitpix = 1;  precision = 'ubit1';
   case   2,
      hdr.image_dimension.bitpix = 8;  precision = 'uint8';
   case   4,
      hdr.image_dimension.bitpix = 16; precision = 'int16';
   case   8,
      hdr.image_dimension.bitpix = 32; precision = 'int32';
   case  16,
      hdr.image_dimension.bitpix = 32; precision = 'float32';
   case  64,
      hdr.image_dimension.bitpix = 64; precision = 'float64';
end;

hdr.image_dimension.dim_un0 = 0;

hdr.image_dimension.pixdim(1) = 0;
hdr.image_dimension.pixdim(2) = pixdim_x;
hdr.image_dimension.pixdim(3) = pixdim_y;
hdr.image_dimension.pixdim(4) = pixdim_z;
hdr.image_dimension.pixdim(5) = 0;
hdr.image_dimension.pixdim(6) = 0;
hdr.image_dimension.pixdim(7) = 0;
hdr.image_dimension.pixdim(8) = 0;


hdr.image_dimension.vox_offset = 0;
hdr.image_dimension.funused1 = 0;
hdr.image_dimension.funused2 = 0;
hdr.image_dimension.funused3 = 0;
hdr.image_dimension.cal_max = 428;
hdr.image_dimension.cal_min = 1;
hdr.image_dimension.compressed = 0;
hdr.image_dimension.verified = 0;
hdr.image_dimension.glmax = max(data(:));
hdr.image_dimension.glmin = min(data(:));

hdr.data_history.descrip = zeros(80, 1, 'int8');
hdr.data_history.aux_file = zeros(24, 1, 'int8');
hdr.data_history.orient = 0;
hdr.data_history.originator = zeros(10, 1, 'int8');
hdr.data_history.generated = zeros(10, 1, 'int8');
hdr.data_history.scannum = zeros(10, 1, 'int8');
hdr.data_history.patient_id = zeros(10, 1, 'int8');
hdr.data_history.exp_date = zeros(10, 1, 'int8');
hdr.data_history.exp_time = zeros(10, 1, 'int8');
hdr.data_history.hist_un0 = zeros(3, 1, 'int8');
hdr.data_history.views = 0;
hdr.data_history.vols_added = 0;
hdr.data_history.start_field = 0;
hdr.data_history.field_skip = 0;
hdr.data_history.omax = 0;
hdr.data_history.omin = 0;
hdr.data_history.smax = 0;
hdr.data_history.smin = 0;


%save the header file
filename = deblank(filename);

fid = fopen(filename, 'w');

fwrite(fid, hdr.header_key.sizeof_hdr, 'int');
fwrite(fid, hdr.header_key.data_type, 'char');
fwrite(fid, hdr.header_key.db_name, 'char');
fwrite(fid, hdr.header_key.extents, 'int');
fwrite(fid, hdr.header_key.session_error, 'short');
fwrite(fid, hdr.header_key.regular, 'char');
fwrite(fid, hdr.header_key.hkey_un0, 'char');

fwrite(fid, hdr.image_dimension.dim, 'short');
fwrite(fid, hdr.image_dimension.unused8, 'short');
fwrite(fid, hdr.image_dimension.unused9, 'short');
fwrite(fid, hdr.image_dimension.unused10, 'short');
fwrite(fid, hdr.image_dimension.unused11, 'short');
fwrite(fid, hdr.image_dimension.unused12, 'short');
fwrite(fid, hdr.image_dimension.unused13, 'short');
fwrite(fid, hdr.image_dimension.unused14, 'short');
fwrite(fid, hdr.image_dimension.datatype, 'short');
fwrite(fid, hdr.image_dimension.bitpix, 'short');
fwrite(fid, hdr.image_dimension.dim_un0, 'short');
fwrite(fid, hdr.image_dimension.pixdim, 'float');
fwrite(fid, hdr.image_dimension.vox_offset, 'float');
fwrite(fid, hdr.image_dimension.funused1, 'float');
fwrite(fid, hdr.image_dimension.funused2, 'float');
fwrite(fid, hdr.image_dimension.funused3, 'float');
fwrite(fid, hdr.image_dimension.cal_max, 'float');
fwrite(fid, hdr.image_dimension.cal_min, 'float');
fwrite(fid, hdr.image_dimension.compressed, 'float');
fwrite(fid, hdr.image_dimension.verified, 'float');
fwrite(fid, hdr.image_dimension.glmax, 'int');
fwrite(fid, hdr.image_dimension.glmin, 'int');

fwrite(fid, hdr.data_history.descrip, 'char');
fwrite(fid, hdr.data_history.aux_file, 'char');
fwrite(fid, hdr.data_history.orient, 'char');
fwrite(fid, hdr.data_history.originator, 'char');
fwrite(fid, hdr.data_history.generated, 'char');
fwrite(fid, hdr.data_history.scannum, 'char');
fwrite(fid, hdr.data_history.patient_id, 'char');
fwrite(fid, hdr.data_history.exp_date, 'char');
fwrite(fid, hdr.data_history.exp_time, 'char');
fwrite(fid, hdr.data_history.hist_un0, 'char');
fwrite(fid, hdr.data_history.views, 'int');
fwrite(fid, hdr.data_history.vols_added, 'int');
fwrite(fid, hdr.data_history.start_field, 'int');
fwrite(fid, hdr.data_history.field_skip, 'int');
fwrite(fid, hdr.data_history.omax, 'int');
fwrite(fid, hdr.data_history.omin, 'int');
fwrite(fid, hdr.data_history.smax, 'int');
fwrite(fid, hdr.data_history.smin, 'int');

fclose(fid);



%%%%save the img file

filename_length = length(filename);
filename(filename_length) = 'g';
filename(filename_length-1) = 'm';
filename(filename_length-2) = 'i';

fid = fopen(filename, 'w');


%     0 None                   (Unknown bit per voxel)   % DT_NONE, DT_UNKNOWN 
%     1 Binary                       (ubit1, bitpix=1)   % DT_BINARY 
%     2 Unsigned char       (uchar or uint8, bitpix=8)   % DT_UINT8, NIFTI_TYPE_UINT8 
%     4 Signed short                 (int16, bitpix=16)  % DT_INT16, NIFTI_TYPE_INT16 
%     8 Signed integer               (int32, bitpix=32)  % DT_INT32, NIFTI_TYPE_INT32 
%    16 Floating point   (single or float32, bitpix=32)  % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
%    32 Complex, 2 float32     (Unsupported, bitpix=64)  % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision (double or float64, bitpix=64)  % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
%   128 uint8 RGB                (Use uint8, bitpix=24)  % DT_RGB24, NIFTI_TYPE_RGB24 
%   256 Signed char          (schar or int8, bitpix=8)   % DT_INT8, NIFTI_TYPE_INT8 
%   511 Single RGB             (Use float32, bitpix=96)  % DT_RGB96, NIFTI_TYPE_RGB96
%   512 Unsigned short              (uint16, bitpix=16)  % DT_UNINT16, NIFTI_TYPE_UNINT16 
%   768 Unsigned integer            (uint32, bitpix=32)  % DT_UNINT32, NIFTI_TYPE_UNINT32 
%  1024 Signed long long             (int64, bitpix=64)  % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long          (uint64, bitpix=64)  % DT_UINT64, NIFTI_TYPE_UINT64 
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
%  1792 Complex128, 2 float64  (Unsupported, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 


switch hdr.image_dimension.dim(1)

	case 4
		if hdr.image_dimension.dim(5) == 1
		
            data_tmp = zeros(hdr.image_dimension.dim(2), hdr.image_dimension.dim(3), hdr.image_dimension.dim(4));
            
			%for i=1:hdr.image_dimension.dim(4)
			%	data_tmp(:,:,i) = squeeze(data(i,:,:))';
			%end;
			
			data_tmp = ipermute(data, [3 2 1]);
			
            %size(data_tmp)
            %hdr.image_dimension.dim(2)*hdr.image_dimension.dim(3)*hdr.image_dimension.dim(3)
            
			data_out = reshape(data_tmp, hdr.image_dimension.dim(2)*hdr.image_dimension.dim(3)*hdr.image_dimension.dim(4), 1);
			clear data_tmp;
			
			fwrite(fid, data_out, precision);
		    
    			clear data_out;
    
		else
      		end
      		
    	case 3
      		% 3-dimensional data
            data_tmp = zeros(hdr.image_dimension.dim(2), hdr.image_dimension.dim(3), hdr.image_dimension.dim(4));
                        
      		%for i=1:hdr.image_dimension.dim(4),
      		%	data_tmp(:,:,i) = squeeze(data(i,:,:))';
      		%end;
      		
      		data_tmp = ipermute(data, [3, 2, 1]);
      		
      		data_out = reshape(data_tmp, hdr.image_dimension.dim(2)*hdr.image_dimension.dim(3)*hdr.image_dimension.dim(4), 1);
      		clear data_tmp;
      		
      		fwrite(fid, data_out, precision);
      		clear data_out;
      		            
    	case 2
      		% 2-dimensional data
end;

fclose(fid);
clear hdr;
clear precision;