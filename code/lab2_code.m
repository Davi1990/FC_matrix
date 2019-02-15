clc
close all
clear

%% Previous steps:
% 1. Coreg Gordon to subject's meanRSFRMI          -> SPM (Hint: use interp: nearest neighbour) 
% 2. Coreg T1_Segmentation to subject's meanRSFRMI -> SPM 


%% Matlab code:

%% Loading Coregistered Data
EPI = '/Volumes/NO NAME/01/fMRI/swaufmri.nii';
parcels = '/Users/ariannamenardi/Documents/MATLAB/bertoldo/atlas/rGordon_parcels_2_T1.nii';
malf = '/Users/ariannamenardi/Documents/MATLAB/bertoldo/atlas/rT1_3D_AnatomicalSegmentation.nii';
parcels_list = '/Users/ariannamenardi/Documents/MATLAB/bertoldo/FC/code/Parcels.xlsx';
moco_params = '/Volumes/NO NAME/01/fMRI/rp_fmri.txt';

EPI     = load_untouch_nii(EPI);
EPI     = double(EPI.img);
malf     = load_untouch_nii(malf);
malf     = double(malf.img);
parcels = load_untouch_nii(parcels);
parcels = double(parcels.img);

%% Creating brain mask (NB : use both segmentations!!!)
mask = (malf + parcels) > 0;

%% Confound Regression:(TAC = extract_ROI_TAC(EPI,parcels,ROI_id))

% Creating confound regression matrix
WM_conf   = extract_ROI_TAC(EPI,malf,[44 45]); % ROI labels: 44 45
CSF_conf  = extract_ROI_TAC(EPI,malf,[51 52]); % ROI labels: 51 52
moco_conf = load(moco_params);                 % load rp*.txt

X = [WM_conf CSF_conf moco_conf];
zX = zscore(X);

% Defining temporal filtering parameters
hp = 0.01; % [Hz]
lp = 0.08; % [Hz]

% Data cleaning
[nR, nC, nS, nVol] = size(EPI);
EPInew = zeros(size(EPI));
for rr = 1: nR
    for cc = 1: nC
        for ss = 1: nS
            if mask(rr,cc,ss)
                voxel = squeeze(EPI(rr,cc,ss,:));
                
                % denoising
                beta  = lscov(zX,voxel);
                noise = zX * beta;
                voxel_cleaned = voxel - noise;
                
                % filtering
                new_voxel = data_filtering(voxel_cleaned,hp,lp);
                
                EPInew(rr,cc,ss,:) = new_voxel;
                
                clear voxel beta noise voxel_cleaned new_voxel
            end
        end
    end
end


%% ROIs time series extraction:
gordonROI = extract_gordon_TACs(EPInew,parcels,parcels_list);

%% Compute FC and zFisher transform
[FC, pval] = corr([gordonROI.TAC]);
zFC = atanh(FC);

%% show FC matrix (according to Gordon well known standard)
view_FC(zFC,parcels_list)
