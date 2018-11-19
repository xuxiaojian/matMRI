%% transform every hdr files to the mat files

clear all
%% don't need to care about the folders 
%% read hdr to mat files
matFolder = '/root/Documents/MATLAB/deepMRI/workSpaceNoNormx2/';
myFolder = '/root/Documents/MATLAB/deepMRI/dataSpaceNoNormx2/Healthy';
cd(myFolder);
d = dir(myFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
totalFolds = length(nameFolds);

% find all the Subj folders and the find the capture and mcnufft folders
% inside
for i = 1:totalFolds
    
    nameSubj = nameFolds(i);
    nameSubj = char(nameSubj);
    cd(strcat(pwd,'/', nameSubj));
    % capture
    csFold = strcat(pwd,'/capture');
    cd (csFold);
    input = dir([pwd, '/*.hdr']);
    hdrFile = input.name;
    matName = strcat('cs', nameSubj);
    genvarname(matName); 
    eval([ matName '= Tools.hdr2mat(hdrFile);']);
    save([matFolder 'cs' nameSubj '.mat'], matName); % second one is the name of the var to save
    cd ..
    % mcnufft
    mcFold = strcat(pwd,'/mcnufft');
    cd (mcFold);
    input = dir(strcat(pwd,'/*.hdr'));
    hdrFile = input.name;
    matName = strcat('mc', nameSubj);
    genvarname(matName); 
    eval([ matName '= Tools.hdr2mat(hdrFile);']); 
    save([matFolder 'mc' nameSubj '.mat'], matName);
    cd ..
    %
    cd ..
    
end

