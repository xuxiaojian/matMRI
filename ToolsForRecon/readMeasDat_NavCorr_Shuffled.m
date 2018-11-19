function [kSpaceData,sampledLineIndices] = readMeasDat_NavCorr_Shuffled(filename,varargin)
% Example:
%
% kSpaceData = readMeasDat('meas.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas.dat',[204 256 70],'Echo',2,'Rep',50);

% fid2 = fopen('headerInfo_WithFlag.txt','w');
ISWrite = 0;
% fid_acq = fopen('acquisitions.txt', 'w');

numberOfPE = 0;
numberOfSlices = 0;
numberOfEchoes = 1;
numberOfRepetitions = 1;

if nargin >= 2
    numberOfPE = varargin{1}(1);
    numberOfSlices = varargin{1}(2);
    if nargin == 4 && strcmpi('echo',varargin{2})
        numberOfEchoes = varargin{3};
    end
    if nargin == 6 && strcmpi('rep',varargin{4})
        numberOfRepetitions = varargin{5};
    end
end

% VB17 RaidFile structure

% FILE *in = fopen("D:\\temp\\meas_orig.dat", "rb");
fid = fopen(filename,'rb');
% int32_t nProtHeaderLen = 0; // read the protocol header length
% fread(&nProtHeaderLen, 1, 4, in);
nProtHeaderLen = double(fread(fid,1,'*int32'));
% fseek(in, nProtHeaderLen - 4, SEEK_CUR);
% fseek(fid,nProtHeaderLen-4,0);
AAA = fread(fid,nProtHeaderLen-4,'uint8');
% Determine the number of lines [#columns and #channels will come from the MDH]
s = '<ParamLong."iNoOfFourierLines">';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
numberOfPE = str2num((neighbor(firstBrace+1:secondBrace-1))); %#ok<*ST2NM>
% Determine the number of partitions
s = '<ParamLong."lPartitions"> ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
numberOfPartitions = str2num((neighbor(firstBrace+1:secondBrace-1)));

amplitudeScaleFactor = 80 * 20 * 131072 / 65536;

amplitudeScaleFactor = amplitudeScaleFactor*20000;

acqEnded = 0;
mdhCounter = 1;
sliceIndices = [];
sampledLineIndices = [];

icount = 0;
while ~acqEnded    
    % MDH Header Information Retrieval Starts Here...
    ulFlagsAndDMALength = fread(fid,1,'*uint32');
    lMeasUID = fread(fid,1,'*int32');
    ulScanCounter = fread(fid,1,'*uint32');
    ulTimeStamp = fread(fid,1,'*uint32');
    ulPMUTimeStamp = fread(fid,1,'*uint32');
    aulEvalInfoMaskMostSig = fread(fid,1,'*uint32');%evaluation info mask field,  first part
    aulEvalInfoMaskLeastSig = fread(fid,1,'*uint32');%evaluation info mask field, second part
    ushSamplesInScan = fread(fid,1,'*uint16'); % # of samples acquired in scan
    ushUsedChannels = fread(fid,1,'uint16');
    ushLine = fread(fid,1,'*uint16');                  % line index                   */
    ushAcquisition = fread(fid,1,'*uint16');           % acquisition index            */
    ushSlice = fread(fid,1,'*uint16');                 % slice index                  */
    ushPartition = fread(fid,1,'*uint16');             % partition index              */
    ushEcho = fread(fid,1,'*uint16');                  % echo index                   */
    ushPhase = fread(fid,1,'*uint16');                 % phase index                  */
    ushRepetition = fread(fid,1,'*uint16');            % measurement repeat index     */
    % fseek(fid,7*2,0); % Other sLoopCounter variables
    ushSet = fread(fid,1,'*uint16');
    ushSeg = fread(fid,1,'*uint16');
    ushIda = fread(fid,1,'*uint16');
    ushIdb = fread(fid,1,'*uint16');
    ushIdc = fread(fid,1,'*uint16');
    ushIdd = fread(fid,1,'*uint16');
    ushIde = fread(fid,1,'*uint16');
    % fseek(fid,4,0); % sCutOff variables
    ushPre =  fread(fid,1,'*uint16');
    ushPost =  fread(fid,1,'*uint16');
    ushKSpaceCenterColumn = fread(fid,1,'*uint16');
    ushCoilSelect = fread(fid,1,'*uint16');
    fReadOutOffCenter = fread(fid,1,'*float32');
    ulTimeSinceLastRF = fread(fid,1,'*uint32');
    ushKSpaceCenterLineNo = fread(fid,1,'*uint16');
    ushKSpaceCenterPartitionNo = fread(fid,1,'*uint16');
    % fseek(fid,8,0); % aushIceProgramPara
    aushIceProgramPara(1) = fread(fid,1,'*uint16');
    aushIceProgramPara(2) = fread(fid,1,'*uint16');
    aushIceProgramPara(3) = fread(fid,1,'*uint16');
    aushIceProgramPara(4) = fread(fid,1,'*uint16');
    % fseek(fid,8,0); % aushFreePara
    aushFreePara(1) = fread(fid,1,'*uint16');
    aushFreePara(2) = fread(fid,1,'*uint16');
    aushFreePara(3) = fread(fid,1,'*uint16');
    aushFreePara(4) = fread(fid,1,'*uint16');
    %fseek(fid,28,0); % sSliceData -> replaced with fSag ... aflQuaternion
    fSag = fread(fid,1,'*float32');
    fCor = fread(fid,1,'*float32');
    fTra = fread(fid,1,'*float32');
    aflQuaternion = fread(fid,4,'*float32');
    %disp([fSag fCor fTra aflQuaternion'])
    ushChannelID = fread(fid,1,'*uint16');
    % fseek(fid,2,0);
    ushPTABPosNeg = fread(fid,1,'*uint16');
    % MDH Header Information Retrieval Ends Here...
    
    mdhBitFields = determineBitFields(aulEvalInfoMaskMostSig);
           
    if mdhCounter==1
        numberOfColumns = double(ushSamplesInScan);
        acquisitionMatrixWidth = numberOfColumns/2;
        evenPhaseCorrectionVector = zeros(numberOfColumns,2);
        numberOfChannels = double(ushUsedChannels);
        evenPhaseCorrectionSoFar = zeros(1,numberOfChannels);
        phaseDifferenceVector = ones(numberOfColumns,numberOfChannels);
        noiseSamples = zeros(numberOfColumns,numberOfChannels);
        if numberOfPE ~=0
            kSpaceData = zeros(numberOfColumns,numberOfPE,numberOfPartitions,numberOfChannels);
        end
%         kSpaceData = [];
    end
    
    % Acquisition ended? [To be checked in the "while" statement above.]
    
    acqEnded = mdhBitFields.MDH_ACQEND;
    
    if acqEnded
        break;
    end
    
    actualData = fread(fid,2*numberOfColumns,'float32=>double');
    
    if ushPartition==0 && ushChannelID==0
        sampledLineIndices = [sampledLineIndices;ushLine];
    end
    
    currentLine = ushLine+1;
    currentPartition = ushPartition + 1;
    currentAcquisition = ushAcquisition+1;
    currentSlice = ushSlice+1;
    currentRepetition = ushRepetition+1;
    currentChannel = ushChannelID+1;
    currentEcho = ushEcho+1;
    currentSegment = ushSeg+1;
    
    %if currentChannel == 1 % && currentLine == 94
    if ~mdhBitFields.MDH_NOISEADJSCAN %&& currentChannel == 1
        if mdhBitFields.MDH_FIRSTSCANINSLICE == 1
            ISWrite = 1;
        end    
%                 disp([fSag fCor fTra aflQuaternion'])
                disp([' Slice: ' num2str(currentSlice) ...
                ' Partition: ' num2str(currentPartition) ...
                ' Line: ' [repmat('0',1,currentLine<10) num2str(currentLine)] ...
                ' Rep: ' num2str(currentRepetition) ...
                ' Seg: ' num2str(currentSegment) ...
                ' Acq: ' num2str(currentAcquisition) ...
                ' Cha: ' num2str(currentChannel) ...
                ' Echo: ' num2str(currentEcho)])
            
    else
        keyboard;
    end
   
    if ISWrite == 1
        if ~mdhBitFields.MDH_PHASCOR
            
            realPart = actualData(1:2:end-1)*amplitudeScaleFactor;
            imaginaryPart = actualData(2:2:end)*amplitudeScaleFactor;
            complexData = realPart+1i*imaginaryPart;
            if mdhBitFields.MDH_REFLECT
                complexData = flipud(complexData);
            else
                complexData = fftshift(fft(    ifft(ifftshift(complexData)).*phaseDifferenceVector(:,currentChannel)    ));
            end            
            
            if ~mdhBitFields.MDH_NOISEADJSCAN                
                try
                    %                     kSpaceData(currentLine,:,currentPartition,currentChannel,currentEcho,currentRepetition) = ...
                    %                     (kSpaceData(currentLine,:,currentPartition,currentChannell,currentEcho,currentRepetition)*(currentAcquisition-1)+...
                    %                     complexData.')/currentAcquisition;
                    kSpaceData(:, currentLine,currentPartition,currentChannel,currentEcho,currentRepetition,currentSlice) = complexData.';
                    %Transposed complex data to a row vector. Furthermore, averaged
                    %kSpaceData over number of acquisitions.
                catch MException %if currentAcquisition = 1, we may not have that portion of
                    %kSpaceData available.
                    %kSpaceData(currentLine,:,currentPartition,currentChannel) = complexData.'; %,currentEcho,currentRepetition
                    kSpaceData(:,currentLine,currentPartition,currentChannel) = complexData.';
                end
            end
        end
    end
    icount = icount + 1;
    mdhCounter = mdhCounter+1;
end


function mdhBitFields = determineBitFields(evalInfo)

bits = num2str(dec2bin(evalInfo));
bits = fliplr([repmat('0',1,32-length(bits)),bits]);
setFlags = find(bits=='1')-1;

mdhBitFields.MDH_ACQEND = any(setFlags==0);
mdhBitFields.MDH_RTFEEDBACK = any(setFlags==1);
mdhBitFields.MDH_HPFEEDBACK = any(setFlags==2);
mdhBitFields.MDH_ONLINE    = any(setFlags==3);
mdhBitFields.MDH_OFFLINE   = any(setFlags==4);
mdhBitFields.MDH_LASTSCANINCONCAT = any(setFlags==8);       % Flag for last scan in concatination
mdhBitFields.MDH_RAWDATACORRECTION = any(setFlags==10);      % Correct the rawadata with the rawdata correction factor
mdhBitFields.MDH_LASTSCANINMEAS = any(setFlags==11);      % Flag for last scan in measurement
mdhBitFields.MDH_SCANSCALEFACTOR = any(setFlags==12);      % Flag for scan specific additional scale factor
mdhBitFields.MDH_2NDHADAMARPULSE = any(setFlags==13);      % 2nd RF exitation of HADAMAR
mdhBitFields.MDH_REFPHASESTABSCAN = any(setFlags==14);      % reference phase stabilization scan
mdhBitFields.MDH_PHASESTABSCAN = any(setFlags==15);      % phase stabilization scan
mdhBitFields.MDH_D3FFT     = any(setFlags==16);      % execute 3D FFT
mdhBitFields.MDH_SIGNREV   = any(setFlags==17);      % sign reversal
mdhBitFields.MDH_PHASEFFT  = any(setFlags==18);      % execute phase fft
mdhBitFields.MDH_SWAPPED   = any(setFlags==19);      % swapped phase/readout direction
mdhBitFields.MDH_POSTSHAREDLINE = any(setFlags==20);      % shared line
mdhBitFields.MDH_PHASCOR   = any(setFlags==21);      % phase correction data
mdhBitFields.MDH_PATREFSCAN = any(setFlags==22);      % additonal scan for PAT reference line/partition
mdhBitFields.MDH_PATREFANDIMASCAN = any(setFlags==23);      % additonal scan for PAT reference line/partition that is also used as image scan
mdhBitFields.MDH_REFLECT   = any(setFlags==24);      % reflect line
mdhBitFields.MDH_NOISEADJSCAN = any(setFlags==25);      % noise adjust scan --> Not used in NUM4
mdhBitFields.MDH_SHARENOW  = any(setFlags==26);      % all lines are acquired from the actual and previous e.g. phases
mdhBitFields.MDH_LASTMEASUREDLINE = any(setFlags==27);      % indicates that the current line is the last measured line of all succeeding e.g. phases
mdhBitFields.MDH_FIRSTSCANINSLICE = any(setFlags==28);      % indicates first scan in slice = any(setFlags==needed for time stamps)
mdhBitFields.MDH_LASTSCANINSLICE = any(setFlags==29);      % indicates  last scan in slice = any(setFlags==needed for time stamps)
mdhBitFields.MDH_TREFFECTIVEBEGIN = any(setFlags==30);      % indicates the begin time stamp for TReff = any(setFlags==triggered measurement)
mdhBitFields.MDH_TREFFECTIVEEND = any(setFlags==31);


function [magnitude,phase] = processSlices(complexSlice,acquisitionMatrixWidth)

magnitude = zeros(acquisitionMatrixWidth,acquisitionMatrixWidth,size(complexSlice,3));
phase = magnitude;
centerWindowFactor = 1/4;

for k = 1:size(complexSlice,3)
    
    complexSlice(:,:,k) = echoShiftCorrection(complexSlice(:,:,k)); % center echo
    complexSlice(:,:,k) = apodization(complexSlice(:,:,k),0.5); % gibb's phenomenon, hanning window to prevent ringing
    fftdata = fftshift(ifft2(ifftshift(complexSlice(:,:,k))));
    fftdata = fftdata(:,acquisitionMatrixWidth/2+1:end-acquisitionMatrixWidth/2);
    magnitude(:,:,k) = abs(fftdata);
    
    complexSlice(:,:,k) = HighPassCmplx(complexSlice(:,:,k), centerWindowFactor);
    fftdata=fftshift(ifft2(ifftshift(complexSlice(:,:,k))));
    fftdata = fftdata(:,acquisitionMatrixWidth/2+1:end-acquisitionMatrixWidth/2);
    phase(:,:,k) = angle(fftdata);
    
end

function saveAnalyzeImage(image,filename)

temp1 = permute(image,[3 1 2]);
% height = size(image,1);
% depth = size(image,3);
analyzeWriter(temp1,1,1,1,16,filename);

%The following function is to high pass filter the phase of the cmplex data by using a window
function filteredKspace = HighPassCmplx(cmplx_ktemp, centerWindowFactor)

[xsize,ysize] = size(cmplx_ktemp);
cmplx_ktemp3=cmplx_ktemp*0;
centerwindowsize = centerWindowFactor*min(xsize,ysize);
halfwsize=centerwindowsize/2;

cmplx_ktemp3(xsize/2-halfwsize:xsize/2+halfwsize-1, ysize/2-halfwsize:ysize/2+halfwsize-1)= ...
    cmplx_ktemp(xsize/2-halfwsize:xsize/2+halfwsize-1, ysize/2-halfwsize:ysize/2+halfwsize-1);

cmplx_itemp2=ifft2(ifftshift(cmplx_ktemp));
cmplx_itemp22=ifft2(ifftshift(cmplx_ktemp3));

cmplx_itemp=cmplx_itemp2.*conj(cmplx_itemp22);
filteredKspace = fftshift(fft2(cmplx_itemp));

function complexSignal = echoShiftCorrection(complexSignal)

%complexSignal = ifftshift(complexSignal);
[maxVal,maxInd] = max(abs(complexSignal(:)));
[xMax,yMax] = ind2sub(size(complexSignal),maxInd);

acquisitionOneSide = size(complexSignal,1);

xshift = acquisitionOneSide/2-xMax+1;
yshift = acquisitionOneSide-yMax+1;

%Echo shift prior to any data analysis to remove global phase shift
complexSignal=circshift(complexSignal,[xshift,yshift]);
%complexSignal = fftshift(complexSignal);

function complexSignal = apodization(complexSignal,skirtingFactor)

% skirtingFactor should be between 0 and 1. The case skirtingFactor = 1
% reduces to the MATLAB hanning() function.

[xsize,ysize] = size(complexSignal);

S_h = ysize*skirtingFactor; % slope length
S_v = xsize*skirtingFactor;

skirt_h = 1/2*(1+cos((1:S_h/2)*2*pi/S_h));
skirt_v = 1/2*(1+cos((1:S_v/2)*2*pi/S_v));

f_h = [fliplr(skirt_h) ones(1,ysize-2*length(skirt_h)+1) skirt_h(1:end-1)];
f_v = [fliplr(skirt_v) ones(1,xsize-2*length(skirt_v)+1) skirt_v(1:end-1)];


complexSignal = repmat(f_v.',1,ysize).*complexSignal;
complexSignal = repmat(f_h,xsize,1).*complexSignal;