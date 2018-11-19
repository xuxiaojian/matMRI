function [projections,maxK,RR,border_largeGA,border_largeGA_filtered_Org] = processCAPTURE_OnServer_VarW_NewQualityMetric_TGV(datFileLocation,firstSlice,lastSlice,nPhases,lambdaFactorWavelet,lambdaFactorInterPhase,lambdaFactorTGV_1stDerivative,bComputeSensitivities,percentW)

bSaveWS = false;

addpath(genpath('/bmr207/nmrgrp/nmr100/ToolsForRecon'))

% Save current directory so that we can go back once all is done.
originalDir = pwd;
%Xiaojian
% cd(originalDir);
% reconDir = './Recons_July18';
% 
% mkdir(reconDir);
 
% Determine folder name and file name
forwardSlashIndices = find(datFileLocation=='/');
folderName = datFileLocation(1:forwardSlashIndices(end)-1);
fileName = datFileLocation(forwardSlashIndices(end)+1:end);

cd(folderName)

%% Read raw data and save as matlab workspaces
wsDir = strrep(fileName,'.dat','');

if bSaveWS % firstSlice==1 && bSaveWS
    % dbquit
    kSpaceData = readMeasDat_NavCorr_Shuffled(fileName);
    mkdir(wsDir)
    cd(wsDir)
    % For Cor and Tra in Cihat's first set, we need to exclude the neck coil
    if ~isempty(strfind(fileName,'MID00434_FID100462_fl3d_vibe_AM_ZeroPhi_Cor_Image')) ...
        || ~isempty(strfind(fileName,'MID00437_FID100465_fl3d_vibe_AM_ZeroPhi_Tra_Image'))
        for k = 2:size(kSpaceData,4)
            eval(['kspace_cha' num2str(k-1) ' = squeeze(kSpaceData(:,:,:,' num2str(k) '));']);
            save(['cha=' num2str(k-1)],['kspace_cha' num2str(k-1)])
            clear(['kspace_cha' num2str(k-1)])
            disp(num2str(k))
        end
    else
        for k = 1:size(kSpaceData,4) % read all the kspace data for the 4th dimension and store them as cha1--k
            eval(['kspace_cha' num2str(k) ' = squeeze(kSpaceData(:,:,:,' num2str(k) '));']);% kspace_chak = = squeeze(kSpaceData(:,:,:,k));
            save(['cha=' num2str(k)],['kspace_cha' num2str(k)]) % chak = kspace_chak
            clear(['kspace_cha' num2str(k)])
            disp(num2str(k))
        end
    end
    clear kSpaceData
else
    if ~isdir(wsDir)
        wsDir = ['meas_' wsDir];
    end
    if ~isdir(wsDir) || length(dir(wsDir))<6 % 6 includes . and .. in addition to the matlab workspaces for at least 4 channels.
        error('Matlab workspaces not ready! Please do a run with firstSlice=1 before running other configurations...');
    else
        cd(wsDir)
    end
end



%% Read protocol from raw data
fid = fopen(['../' fileName]);
nProtHeaderLen = double(fread(fid,1,'*int32'));
AAA = fread(fid,nProtHeaderLen-4,'uint8');

% Read center partition no from first MDH
fseek(fid,76,0);
ushKSpaceCenterLineNo = fread(fid,1,'*uint16'); %#ok<NASGU>
ushKSpaceCenterPartitionNo = fread(fid,1,'*uint16');
fclose(fid);

s = 'sCoilElementID.tElement';
ind = strfind(AAA',s);
coilIDs = cellstr(char(AAA(repmat(ind',1,3)+length(s)+repmat(4:6,length(ind),1))));
numberOfChannels = length(coilIDs)/2; % ASCCONV appears twice in the file header...
coilIDs = coilIDs(1:numberOfChannels);
sortedSPCoils = sort(coilIDs(strncmp(coilIDs,'SP',2)));
selectedSPCoil = sortedSPCoils(floor((length(sortedSPCoils)+1)/2));
middleSPIndex = find(strcmp(coilIDs,selectedSPCoil{1}));

s = '<ParamLong."lRadialViews"> ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
ntviews = str2num((neighbor(firstBrace+1:secondBrace-1)));
ntviews = min(ntviews,1000); % Use first 2000 if #Views>2000

if any(strfind(folderName,'Subj3/Cor'))
    % Neck coil was excluded in Cihat's data. So, readjust:
    coilIDs = coilIDs(2:end);
    middleSPIndex = find(strcmp(coilIDs,'SP2'));
    numberOfChannels = length(coilIDs);
end

if any(strfind(folderName,'Subj3/Tra')) % Exception added on May 16, 2017
    % Neck coil was excluded in Cihat's data. So, readjust:
    coilIDs = coilIDs(2:end);
    middleSPIndex = find(strcmp(coilIDs,'SP2'));
    numberOfChannels = length(coilIDs);
end

if any(strfind(folderName,'Subj10/Sag'))
    % Exclude SP4 in Subj10's Sag data because it is very low, leading to
    % very small signal and hence corrupting normalization. Scanner did not
    % include it in Cor or Tra.
    coilIDs(strcmp(coilIDs,'SP4')) = [];
    middleSPIndex = find(strcmp(coilIDs,'SP2'));
    numberOfChannels = length(coilIDs);
end

a = dir('fl3d_vibe_radial_sag_cha=*.mat');
for k = 1:length(a)
    movefile(a(k).name,strrep(a(k).name,'fl3d_vibe_radial_sag_',''));
end

% Use the first 2000 lines
% [make sure nPhases>=10]
load cha=1.mat
K(:,:,1) = kspace_cha1(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1); %#ok<*NODEF>
clear kspace_cha1;
disp(['Done with Channel ' num2str(1)])
load cha=2.mat
K(:,:,2) = kspace_cha2(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
clear kspace_cha2;
disp(['Done with Channel ' num2str(2)])
load cha=3.mat
K(:,:,3) = kspace_cha3(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
clear kspace_cha3;
disp(['Done with Channel ' num2str(3)])
load cha=4.mat
K(:,:,4) = kspace_cha4(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
clear kspace_cha4;
disp(['Done with Channel ' num2str(4)])
if numberOfChannels == 5
    load cha=5.mat
    K(:,:,5) = kspace_cha5(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
    clear kspace_cha5;
    disp(['Done with Channel ' num2str(5)])
end
if numberOfChannels == 6
    load cha=5.mat
    K(:,:,5) = kspace_cha5(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
    clear kspace_cha5;
    disp(['Done with Channel ' num2str(5)])
    
    load cha=6.mat
    K(:,:,6) = kspace_cha6(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
    clear kspace_cha6;
    disp(['Done with Channel ' num2str(6)])
end
if numberOfChannels == 7
    load cha=5.mat
    K(:,:,5) = kspace_cha5(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
    clear kspace_cha5;
    disp(['Done with Channel ' num2str(5)])
    
    load cha=6.mat
    K(:,:,6) = kspace_cha6(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
    clear kspace_cha6;
    disp(['Done with Channel ' num2str(6)])
    
    load cha=7.mat
    K(:,:,7) = kspace_cha7(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:ntviews,1);
    clear kspace_cha7;
    disp(['Done with Channel ' num2str(7)])
end

% W = repmat([zeros(size(K,1)/4,1);hamming(size(K,1)/2);zeros(size(K,1)/4,1)],1,size(K,2));
% W = repmat([zeros(size(K,1)*3/8,1);hamming(size(K,1)/4);zeros(size(K,1)*3/8,1)],1,size(K,2));

nonzeroWidth = size(K,1)*percentW/100;
zeroWidth = size(K,1)-nonzeroWidth;

W = repmat([zeros(zeroWidth/2,1);hamming(nonzeroWidth);zeros(zeroWidth/2,1)],1,size(K,2));

% New quality metric needs the following line commented.
% projections = fftshift(ifft(ifftshift(W.*K(:,:,middleSPIndex),1),[],1),1); 

% Read number of measured partitions from raw data
s = '<ParamLong."NParMeas">';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
nParMeas = str2num((neighbor(firstBrace+1:secondBrace-1))); %#ok<*ST2NM>

% Determine the number of slices to be reconstructed
s = '<ParamLong."lImagesPerSlab"> ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
numberOfSlices = str2num((neighbor(firstBrace+1:secondBrace-1)));

if lastSlice>numberOfSlices
    lastSlice = numberOfSlices;
    if firstSlice>lastSlice
        firstSlice=lastSlice;
    end
end

% Determine the FOV in mm
s = '<ParamDouble."RoFOV">  { <Precision> 16 ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
secondBrace = find(neighbor=='}',1);
FOV = str2num(neighbor(1:secondBrace-1));

% Determine the number of points per spoke
s = '<ParamLong."NImageLins"> ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
ndata =   2   *  str2num((neighbor(firstBrace+1:secondBrace-1)));

% Determine the slice thickness
s = '<ParamArray."SliceThickness"> ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+250)');
openingBraces = find(neighbor=='{');
closingBraces = find(neighbor=='}');
sliceThickness = (1/numberOfSlices)*str2num((neighbor(openingBraces(3)+1:closingBraces(2)-1)));

% Determine the orientation
s = '<ParamLong."lSag"> ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
isSAG = ~isempty(str2num((neighbor(firstBrace+1:secondBrace-1))));

s = '<ParamLong."lCor"> ';
ind = strfind(AAA',s);
neighbor = char(AAA(ind(1)+length(s):ind(1)+length(s)+50)');
firstBrace = find(neighbor=='{');
secondBrace = find(neighbor=='}');
isCOR = ~isempty(str2num((neighbor(firstBrace+1:secondBrace-1))));

isTRA = ~isSAG & ~isCOR;

% Read TR from raw data
indTR = strfind(AAA','alTR[0]');
indEquals = strfind(AAA(indTR(1):indTR(1)+100)','=');
indContrasts = strfind(AAA(indTR(1):indTR(1)+100)','lContrasts');
TR = str2double(char(AAA(indTR(1)+indEquals(1):indTR(1)+indContrasts-2)'))/1000; 

tic;
T = TR*nParMeas*1e-3+19e-3; % 19e-3 is for Q-fat sat
Fs = 1/T;

% N = size(projections,2);
N = size(K,2);
f = -0.5*Fs:Fs/N:0.5*Fs-Fs/N;

% RR = zeros(100,1);
% for k = 1:100
%     [~,border_largeGA] = max(real(projections*exp(-1j*2*pi*k/100)),[],1);
%     magSpectrum = abs(fftshift(fft(border_largeGA-mean(border_largeGA))));
%     RR(k) = sum(magSpectrum(abs(f)<0.5 & abs(f)>0.1))./sum(magSpectrum(abs(f)>0.8));
% end
% [~,maxK] = max(RR);
% k = maxK;

% New quality metric block begin
RR = zeros(100,numberOfChannels);
for ci = 1:numberOfChannels
    projections = fftshift(ifft(ifftshift(W.*K(:,:,ci),1),[],1),1);    
    for k = 1:100
        [~,border_largeGA] = max(real(projections*exp(-1j*2*pi*k/100)),[],1);
        magSpectrum = abs(fftshift(fft(border_largeGA-mean(border_largeGA))));
        RR(k,ci) = sum(magSpectrum(abs(f)<0.5 & abs(f)>0.1))./sum(magSpectrum(abs(f)>0.8))/range(border_largeGA);
    end
end
[~,maxInd] = max(RR(:));
[maxK,ci] = ind2sub(size(RR),maxInd);
disp('List of channels:')
disp(char(coilIDs))
disp(['Selected channel: ' coilIDs{ci}])
projections = fftshift(ifft(ifftshift(W.*K(:,:,ci),1),[],1),1);
k = maxK;
% New quality metric block end

[~,border_largeGA] = max(real(projections*exp(-1j*2*pi*k/100)),[],1);

b = fir1(11,1/(Fs/2));
a = 1;
border_largeGA_filtered = filtfilt(b,a,border_largeGA);

border_largeGA_filtered_SG = filtfilt(b,a,sgolayfilt(border_largeGA,1,5));

border_largeGA_filtered([1:11 end-10:end]) = border_largeGA_filtered_SG([1:11 end-10:end]);

border_largeGA_filtered_Org = border_largeGA_filtered;

%% Outlier removal

indicesToFurtherDiscard = [];

border_largeGA_filtered(indicesToFurtherDiscard) = [];

[~,sortingInd] = sort(border_largeGA_filtered);

ntviews2 = length(sortingInd); 
nspokes = length(sortingInd)/nPhases; 
nt = floor(ntviews2/nspokes);
binningResult = zeros(length(sortingInd),1);
for ind_group = 0 : nt - 1
    spokes_ind = sortingInd(nspokes * ind_group + 1 : nspokes * ind_group + nspokes);
    binningResult(spokes_ind) = ind_group+1;
end

fidTable = fopen('CAPTURE_Table.txt','w');
fprintf(fidTable,'nPhases = %d\n',nPhases);
fprintf(fidTable,'ntviews = %d\n',ntviews);
fprintf(fidTable,'nParMeas = %d\n',nParMeas);
fprintf(fidTable,'TR = %1.4f ms\n',TR);
fprintf(fidTable,'T_FatSat = %1.4f ms\n\n',19);
for k = 1:max(max(nPhases,10),nPhases*ceil(10/nPhases))
    fprintf(fidTable,'Stack %d --> Discarded\n',k);
end
for k = 1:length(binningResult)
    fprintf(fidTable,'%d %1.4f\n',binningResult(k),border_largeGA_filtered(k));
end
fclose(fidTable);

binningTime = toc;
disp(['Binning done in ' num2str(binningTime) ' seconds'])

% save(['../HammingWS_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_np4W=' num2str(percentW*640/100) '_' num2str(nPhases) 'Phs_NQM.mat'],'projections','maxK','ci','coilIDs','RR','border_largeGA','border_largeGA_filtered_Org','indicesToFurtherDiscard');

%%

ncha = numberOfChannels;
for icha = 1 : ncha        
    load(['cha=', num2str(icha),'.mat']);  
    disp(['Loaded data for Channel ' num2str(icha)])
end

cd ..

%the number of slices after CS reconstruction
nslc_f = numberOfSlices;% 96
%the number of slices acquired/before CS reconstruction
nslc = nParMeas; % 38;

kSpaceData = zeros(ndata, ntviews, nslc, ncha);

kSpaceData(:, :, :, 1) = kspace_cha1(:,1:ntviews,1:nslc);
clear kspace_cha1
kSpaceData(:, :, :, 2) = kspace_cha2(:,1:ntviews,1:nslc);
clear kspace_cha2
kSpaceData(:, :, :, 3) = kspace_cha3(:,1:ntviews,1:nslc);
clear kspace_cha3
kSpaceData(:, :, :, 4) = kspace_cha4(:,1:ntviews,1:nslc);
clear kspace_cha4
if numberOfChannels == 5
    kSpaceData(:, :, :, 5) = kspace_cha5(:,1:ntviews,1:nslc);
    clear kspace_cha5
end
if numberOfChannels == 6
    kSpaceData(:, :, :, 5) = kspace_cha5(:,1:ntviews,1:nslc);
    clear kspace_cha5
    kSpaceData(:, :, :, 6) = kspace_cha6(:,1:ntviews,1:nslc);
    clear kspace_cha6
end
if numberOfChannels == 7
    kSpaceData(:, :, :, 5) = kspace_cha5(:,1:ntviews,1:nslc);
    clear kspace_cha5
    kSpaceData(:, :, :, 6) = kspace_cha6(:,1:ntviews,1:nslc);
    clear kspace_cha6
    kSpaceData(:, :, :, 7) = kspace_cha7(:,1:ntviews,1:nslc);
    clear kspace_cha7
end

% the slice with kz=0
% islc_center = 15; % 1-based. In the MDH, the central partition no is 14.
islc_center = double(ushKSpaceCenterPartitionNo)+1;

% IMPORTANT!!! Speed up processing by allocating only needed memory.
% Deleting chunks of data from the 4D matrix takes a long time!
selectedSpokes = 1:ntviews;
selectedSpokes(1:max(max(nPhases,10),nPhases*ceil(10/nPhases))) = []; % Wait for stead-state [make sure nPhases>=10]
selectedSpokes(indicesToFurtherDiscard) = [];

kSpaceData_t = zeros(ndata, length(selectedSpokes), nslc_f, ncha);
maxAbsKspaceData_t = 0;

kdata_f = zeros(nslc_f, ndata);

diff_kdata = (nslc_f / 2  - islc_center)+1;
for icha = 1 : ncha
    disp(['iFFTz for Channel ' num2str(icha)])
    kdata = squeeze(kSpaceData(:, :, :, icha));
    for ispoke = 1 : length(selectedSpokes)
        % disp(['Spoke ' num2str(selectedSpokes(ispoke))])
        % --- SKIP NAVIGATOR ---
        kdata_f(1 : diff_kdata+1, :) = 0; % --> +1 for navigator
        for islc = 2 : nslc % Started from 2 due to navigator
            kdata_f(islc + diff_kdata, :) = kdata(:,selectedSpokes(ispoke),islc);
        end    
        kdata_f((nslc + diff_kdata + 1) : nslc_f, :) = 0;
        %one-dimensional fft
        kdata_f = flip(fftshift(ifft(ifftshift(kdata_f,1), [], 1),1),1);
        for islc = 1 : nslc_f
            kSpaceData_t(:, ispoke, islc, icha) = kdata_f(islc, :);
        end
        if max(abs(kdata_f(:)))>maxAbsKspaceData_t
            maxAbsKspaceData_t = max(abs(kdata_f(:)));
        end
    end
end
clear kdata_f kdata kSpaceData

%% 
%calculate the density function for nufft reconstruction of all spokes
disp('Getting radial spokes...')
K = GenerateRadialSpokesNew(ntviews, ndata,1:ntviews);
W = abs(-1:2/ndata:1-2/ndata)';
w = zeros(ndata, ntviews);
for iw = 1 : ntviews
    w(:, iw) = W(:);
end
N = [ndata, ndata];
[xx,yy] = meshgrid(-1:2/N(1):1-2/N(1));
ph = 1;
disp('Getting FT operator...')

% Throw away first 10 samples
K = K(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:end);
w = w(:,max(max(nPhases,10),nPhases*ceil(10/nPhases))+1:end);

% Also, discard the indices that fall outside the 1st and 99th percentiles
K(:,indicesToFurtherDiscard)=[];
w(:,indicesToFurtherDiscard)=[];

%%
FT = NUFFT(K,1, ph, 0,N, 2);
tmp=zeros(N);
tmp(end/2+1,end/2+1)=1;
tmp=FT'*(w.*(FT*tmp));
%the density function
w = w/max(abs(tmp(:)));

kSpaceData_t = kSpaceData_t/maxAbsKspaceData_t; % TGV

%% TGV Parameters
if bComputeSensitivities
%     % Xiaojian
%     cd(reconDir);
    
    senseIter = 100;        % usually there is no need to change this

    if ~isdir(['MCoe640_WS_TGV_' fileName([6:24 end-12:end-4]) '_NQM'])
        mkdir(['MCoe640_WS_TGV_' fileName([6:24 end-12:end-4]) '_NQM']);
    end

    imgNUFFT = zeros(ndata / 2, ndata / 2, lastSlice-firstSlice+1);

    for slcCounter = 1:lastSlice-firstSlice+1 

        sliceIndex = firstSlice+slcCounter-1;

        imgSens = zeros(ndata, ndata, ncha); % TGV

        disp(['Sensitivity estimation for slice ' num2str(sliceIndex)])    

        img_dc = zeros(ndata / 2, ndata / 2, ncha);

        for icha = 1 : ncha

            data = squeeze(kSpaceData_t(:, :, sliceIndex,icha));

            %do nufft reconstruction
            im_dc = FT'*(data.*w);

            %take the 320x320 points from a 640x640 image
            img_dc(:, :, icha) = im_dc(ndata/2-ndata/4:ndata/2+ndata/4 - 1, ndata/2-ndata/4:ndata/2+ndata/4 - 1);
            imgSens(:,:,icha) = h1_l2_2D_pd(im_dc, sqrt(w).*kSpaceData_t(:, :, sliceIndex,icha), FT, w, 1000, senseIter, 0.001); % TGV

            % smooth sensitivity estimate [TGV]
            mask = [0 1 0; 1 4 1; 0 1 0]/8;
            for j=1:50
                imgSens(:,:,icha) = conv2(imgSens(:,:,icha),mask,'same');
            end

%             % Cihat's correction begins here...
%             imgSens(:,:,icha) = imgSens(:,:,icha)/max(max(abs(imgSens(:,:,icha)))).^2;
%             % Cihat's correction ends here...
        end

        % TGV
        imgSensSOS = sqrt(sum(abs(imgSens).^2,3));
        imgSens = imgSens./repmat(imgSensSOS,[1 1 ncha]);
        imgSens(isnan(imgSens))=0;

        imgNUFFT(:,:,slcCounter) = sqrt(sum(abs(img_dc).^2,3));
       
        
        MCoe = imgSens;
        % If no MCoe data, then here should be '/MCoe_S' instead of _NQM/MCoe_S';
        save(['MCoe640_WS_TGV_' fileName([6:24 end-12:end-4]) '_NQM/MCoe_S' num2str(sliceIndex) '.mat'],'MCoe')
    end
    
    % Xiaojian comment

%     if isCOR
%         actualToAnalyze(flip(flip(imgNUFFT,1),2),['Cor_Unbinned_NUFFT_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_NQM_T=400.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
%     elseif isSAG
%         actualToAnalyze(flip(imgNUFFT,2),['Sag_Unbinned_NUFFT_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_NQM_T=400.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
%     elseif isTRA
%         actualToAnalyze(flip(permute(imgNUFFT,[2 1 3]),1),['Tra_Unbinned_NUFFT_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_NQM_T=400.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
%     end
end


%% CAUTION: ntviews2 is further reduced after the percentiles operation
ntviews2 = length(sortingInd); 
nspokes = length(sortingInd)/nPhases; 
ndata = 640; 
nt = floor(ntviews2/nspokes);

K_traj = zeros(ndata, ntviews2);
w_traj = zeros(ndata, ntviews2);
for ind_group = 0 : nt - 1        
    spokes_ind = sortingInd(nspokes * ind_group + 1 : nspokes * ind_group + nspokes);            
    K2 = K(:,spokes_ind);
    K2 = K2/max(abs(K2(:)))/2;    
    K_traj(:, ind_group * nspokes + 1: (ind_group + 1) * nspokes) = K2;
    disp('Calculating density compensation...')
    
    % Cihat's correction begins here
%     w = voronoidens(K2(:)); % calculate voronoi density compensation function
%     w = w/max(w(:));
%     w = reshape(w,ndata,nspokes);
     
    [CC,IAA,ICC] = unique(K2(:));  
    w = voronoidens(CC);  
    % Duplicate density for previously-removed points [i.e. DC points]
    w = w(ICC);     
    w = reshape(w,ndata,nspokes);
    % Finally, distribute the weight of DC to all DC points equally.
    w(ndata/2+1,:) = w(ndata/2+1,:)/nspokes;
%     % Normalize as before ---> THIS NORMALIZATION IS WRONG - IT CAUSES INTER-PHASE INTENSITY VARIATIONS
%     w = w/max(w(:));
    % Cihat's correction ends here

    w_traj(:, ind_group * nspokes + 1: (ind_group + 1) * nspokes) = w;      
    disp(['Done with ind_group = ' num2str(ind_group)])
end

%%
w = w_traj;

recon_MCNUFFT = zeros(ndata/2,ndata/2,nPhases,lastSlice-firstSlice+1);
recon_CS = zeros(ndata/2,ndata/2,nPhases,lastSlice-firstSlice+1);
    
for slcCounter = 1:lastSlice-firstSlice+1
    
    sliceIndex = firstSlice+slcCounter-1;
    
    disp(['Processing slice ' num2str(sliceIndex)])
% Uncomment by Xiaojian
    load(['MCoe640_WS_TGV_' fileName([6:24 end-12:end-4]) '_NQM/MCoe_S' num2str(sliceIndex) '.mat'])
%     load(['MCoe640_WS_TGV_' fileName([6:24 end-12:end-4]) '/MCoe_S' num2str(sliceIndex) '.mat'])
    
    b1 = MCoe; 
    b1 = b1/max(abs(b1(:)));
    kdata = squeeze(kSpaceData_t(:,:,sliceIndex,:));
    [nx,~,nc]=size(kdata);
    
    [~,invSortingInd] = sort(sortingInd);
    
    for ch=1:nc
        kdata(:,:,ch)=kdata(:,:,ch).*sqrt(w(:,invSortingInd));
    end        
    k = K_traj;

    % sort the data into a time-series
    kdatau = zeros(nx,nspokes,nc,nt);
    ku = zeros(nx,nspokes,nt);
    wu = zeros(nx,nspokes,nt);
    for ii=1:nt
        kdatau(:,:,:,ii)=kdata(:,sortingInd((ii-1)*nspokes+1:ii*nspokes),:);
        ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
        wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
    end
    % multicoil NUFFT operator
    disp('Performing multicoil NUFFT')
    param.E=MCNUFFT(ku,wu,b1);
    % undersampled data
    param.y=kdatau;
    clear kdata kdatau k ku wu %w
    
    % nufft recon
    recon_nufft=param.E'*param.y;
    recon_MCNUFFT(:,:,:,slcCounter) = recon_nufft(ndata/4+1:3*ndata/4,ndata/4+1:3*ndata/4,:);
    
    
    % xiaojian uncomment on Nov/19/2018
    % can be commet below
        param.W = Wavelet_CCIR('Daubechies',4,4);	% Wavelet
        param2.W = TV_Temp(); % CompositeL1    
        param3.W = TGV_dx_CCIR(); % TGV
        param4.W = TGV_dy_CCIR(); % TGV
        param5.W = TGV_d2x_CCIR(); % TGV
        param6.W = TGV_d2y_CCIR(); % TGV

        param.lambda=lambdaFactorWavelet*max(abs(recon_nufft(:)));
        param2.lambda=lambdaFactorInterPhase*max(abs(recon_nufft(:)));
        param3.lambda=lambdaFactorTGV_1stDerivative*max(abs(recon_nufft(:)));
        param4.lambda=lambdaFactorTGV_1stDerivative*max(abs(recon_nufft(:)));

        lambdaFactorTGV_2nd = 2*lambdaFactorTGV_1stDerivative; % Based on Knoll's paper
        param5.lambda=lambdaFactorTGV_2nd*max(abs(recon_nufft(:)));
        param6.lambda=lambdaFactorTGV_2nd*max(abs(recon_nufft(:)));

        param.nite = 8;
        param.display=1;
        fprintf('\n GRASP reconstruction \n')    
        recon_cs=recon_nufft;
        for n=1:2
            recon_cs = CSL1NlCg_CompositeL1_TGV(recon_cs,param,param2,param3,param4,param5,param6);%
        end
        recon_CS(:,:,:,slcCounter) = recon_cs(ndata/4+1:3*ndata/4,ndata/4+1:3*ndata/4,:);%
        disp(['Done with Slice ' num2str(sliceIndex) ' for ' strrep(folderName,'/bmr207/nmrgrp/nmr100/','')])
    
    % can be comment above 
end


%% xiaojian
% cd(reconDir);

if isCOR
    actualToAnalyze(flip(flip(abs(squeeze(recon_MCNUFFT))*1e6,1),2),['Cor_Binned_MCNUFFT_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_np4W=' num2str(percentW*640/100) '_' num2str(nPhases) 'Phs_NQM_' coilIDs{ci} '_TGV_T=1000_NoNorm.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
elseif isSAG
    actualToAnalyze(flip(abs(squeeze(recon_MCNUFFT)),2)*1e6,['Sag_Binned_MCNUFFT_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_np4W=' num2str(percentW*640/100) '_' num2str(nPhases) 'Phs_NQM_' coilIDs{ci} '_TGV_T=1000_NoNorm.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
elseif isTRA
    actualToAnalyze(flip(permute(abs(squeeze(recon_MCNUFFT))*1e6,[2 1 3 4]),1),['Tra_Binned_MCNUFFT_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_np4W=' num2str(percentW*640/100) '_' num2str(nPhases) 'Phs_NQM_' coilIDs{ci} '_TGV_T=1000_NoNorm.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
end

if isCOR
    actualToAnalyze(flip(flip(abs(squeeze(recon_CS))*1e6,1),2),['Cor_Binned_CS_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_np4W=' num2str(percentW*640/100) '_' num2str(nPhases) 'Phs_NQM_' coilIDs{ci} '_TGV_Lm=' strrep(num2str(lambdaFactorInterPhase),'.','p') '_Ls=' strrep(num2str(lambdaFactorTGV_1stDerivative),'.','p') '_T=400.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
elseif isSAG
    actualToAnalyze(flip(abs(squeeze(recon_CS)),2)*1e6,['Sag_Binned_CS_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_np4W=' num2str(percentW*640/100) '_' num2str(nPhases) 'Phs_NQM_' coilIDs{ci} '_TGV_Lm=' strrep(num2str(lambdaFactorInterPhase),'.','p') '_Ls=' strrep(num2str(lambdaFactorTGV_1stDerivative),'.','p') '_T=400.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
elseif isTRA
    actualToAnalyze(flip(permute(abs(squeeze(recon_CS))*1e6,[2 1 3 4]),1),['Tra_Binned_CS_' fileName([6:24 end-12:end-4]) '_S' num2str(firstSlice) '-' num2str(lastSlice) '_np4W=' num2str(percentW*640/100) '_' num2str(nPhases) 'Phs_NQM_' coilIDs{ci} '_TGV_Lm=' strrep(num2str(lambdaFactorInterPhase),'.','p') '_Ls=' strrep(num2str(lambdaFactorTGV_1stDerivative),'.','p') '_T=400.hdr'],FOV/(ndata/2),FOV/(ndata/2),sliceThickness);
end

cd(originalDir)
