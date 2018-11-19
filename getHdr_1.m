clear all

% doData
fileName01 = 'Subj01/meas_MID00384_FID100412_fl3d_vibe_AM_ZeroPhi_Sag_Image.dat';
fileName02 = 'Subj02/meas_MID00405_FID100433_fl3d_vibe_AM_ZeroPhi_Sag_Image.dat';
fileName03 = 'Subj03/meas_MID00431_FID100459_fl3d_vibe_AM_ZeroPhi_Sag_Image.dat';
fileName04 = 'Subj04/meas_MID00631_FID100655_fl3d_vibe_AM_ZeroPhi_Sag_Image.dat';
fileName05 = 'Subj05/meas_MID00081_FID100745_fl3d_vibe_AM_ZeroPhi_Sag_Image.dat';
fileName06 = 'Subj06/meas_MID00118_FID105583_CAPTURE_GA_Sag_Image.dat';
fileName07 = 'Subj07/meas_MID00197_FID106044_CAPTURE_GA_Sag_Image.dat';
fileName08 = 'Subj08/meas_MID00220_FID106067_CAPTURE_GA_Sag_Image.dat';
fileName09 = 'Subj09/meas_MID00655_FID106502_CAPTURE_GA_Sag_Image.dat';
fileName10 = 'Subj10/meas_MID00301_FID114358_CAPTURE_GA_Sag_Image.dat';
% changed all the Mcoefiles; add Subj3,9,10
%% need to be changed every time!!!!!
fileName = fileName10; %% change this every time

%%
pathName = '/root/Documents/MATLAB/deepMRI/dataSpaceNoNormx2/Healthy/';
lastSlice = 96; % check this every time;
phase = 10; % check this every time;
fullName = [pathName,fileName];
processCAPTURE_OnServer_VarW_NewQualityMetric_TGV(fullName,1,lastSlice,phase,0,0.1,0.00125,false,12.5);





















