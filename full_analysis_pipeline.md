This is the analysis code for "Transdiagnostic neural variability deficits in three major psychiatric disorders"
The preprocessing of the T1w/rest-fMRI is not shown here, which can be easily done with publicly avilable tool like fMRIprep/nilearn/DPABI
Main contents include calculating SDBOLD, performing 4-group partial least squares correlation (PLSC) (https://github.com/MIPLabCH/myPLS) 
then calculating the Bootstraping ratio (BSR) map, and decode the map with neuromaps (https://github.com/netneurolab/neuromaps)

## Step1: calculating the SDBOLD for each individual
This is depending on MATLAB with spm and freesurfer installed

```
%% initializing
data_dir = './fMRI_surf/';
lh_data_list = dir([data_dir 'fsaverage4/lh*']);
rh_data_list = dir([data_dir 'fsaverage4/rh*']);
cd([data_dir 'fsaverage4']);
surf_lh=['./lh.pial.fsaverage4.surf.gii'];
surf_rh=['./rh.pial.fsaverage4.surf.gii'];

% get mask
MaskData_lh_gheader=gifti(./'lh.cortex.fsaverage4.func.gii');
MaskData_rh_gheader=gifti(./'rh.cortex.fsaverage4.func.gii');
MaskData_lh=double(logical(MaskData_lh_gheader.cdata)); MaskDataOneDim_lh=reshape(MaskData_lh,1,[]);
MaskData_rh=double(logical(MaskData_rh_gheader.cdata)); MaskDataOneDim_rh=reshape(MaskData_rh,1,[]);

%% read data and do calculating

for idx = 1:length(lh_data_list)
    
    %% our fMRI data in fsaverage space is in nifti format, so the dimensions need to be reshaped, the
    %% fMRI data in gifit format (e.g. *.func.gii) shoud be directly read as a 2-D matrix by 'gifti' function
    %% each row is a vertex and each coloum is a time point 
    
    mri_lh = MRIread(lh_data_list(idx).name);
    mri_rh = MRIread(rh_data_list(idx).name);
    data_lh = reshape(mri_lh.vol,[size(mri_lh.vol,1)*size(mri_lh.vol,2)*size(mri_lh.vol,3), ...
        size(mri_lh.vol,4)]);
    data_rh = reshape(mri_rh.vol,[size(mri_rh.vol,1)*size(mri_rh.vol,2)*size(mri_rh.vol,3), ...
        size(mri_rh.vol,4)]);
    
    % calculate normalized SD bold
    sd_lh(:, idx) = (std(data_lh'))';
    sd_rh(:, idx) = (std(data_rh'))';
    sd_all_temp = [sd_lh(MaskData_lh==1, idx); sd_rh(MaskData_rh==1, idx)];
    mean_sd = mean(sd_all_temp);
    std_sd = std(sd_all_temp);
    
    sd_lh_z(:, idx) = (sd_lh(:, idx)-mean_sd)./std_sd;
    sd_rh_z(:, idx) = (sd_rh(:, idx)-mean_sd)./std_sd;
    
    % we dont want the medial wall vertex
    sd_lh_z(MaskData_lh==0, idx) = 0;
    sd_rh_z(MaskData_rh==0, idx) = 0;
    clear sd_all_temp mean_sd std_sd
end

sd_data_z = [sd_lh_z, sd_rh_z]
```

## Step 2: perform the PLSC with myPLS toolbox (https://github.com/MIPLabCH/myPLS)
the original myPLS toolbox dont support group number higher than 3, so we need to generate the contrast matrix
as the Y (behavior) input and use the "behavior" type PLSC

```
% generate contrast
% here 0 for controls, 1 for bipolar manic, 2 for bipolar depression, 3 for schizophrenia, 4 for depression

contrast=zeros(498,4);

contrast(diagnosis~=0,1)=1./sum(diagnosis~=0);
contrast(diagnosis==0,1)=-1./sum(diagnosis==0);
contrast(:,1)=orth(contrast(:,1));

contrast(diagnosis==2,2)=1;
contrast(diagnosis==1,2)=1;
contrast(diagnosis==4,2)=1;
contrast(contrast(:,2)==1,2)=1/sum(contrast(:,2)==1);
contrast(diagnosis==3,2)=-1./sum(diagnosis==3);
contrast(:,2)=orth(contrast(:,2));

contrast(diagnosis==4,3)=1;
contrast(diagnosis==2,3)=1;
contrast(contrast(:,3)==1,3)=1/sum(contrast(:,3)==1);
contrast(diagnosis==1,3)=-1./sum(diagnosis==1);
contrast(:,3)=orth(contrast(:,3));

contrast(diagnosis==2,4)=-1/sum(diagnosis==2);
contrast(diagnosis==4,4)=1/sum(diagnosis==4);
contrast(:,4)=orth(contrast(:,4));

% load data of age, gender and EDUY
% a vector contain one type of data
load demographic_data

% regress effect of demopraphic data

Y = sd_data_z; % generate from previous step 1
NV_reg = [age, gender, EDUY];
Y_CN   = Y(diagnosis==0,:);
X_CN   = [ones(size(Y_CN,1),1) NV_reg(diagnosis==0,:)];
mean_CN = mean(Y_CN,1);
% add a diagonal matrix with small values to prevent singular matrix
% problem
b_CN = (X_CN'*X_CN + eye(size(X_CN,2))*1e-8)\(X_CN'*Y_CN);
% Regress out regressors for all subjects by computing Y - X*b_CN +
% mean(Y_CN) for each unique ROI-ROI pair
X = [ones(size(Y,1),1) NV_reg];
sd_data_z_reg  = bsxfun(@plus, Y-X*b_CN, mean_CN);
clear NV_reg X X_CN Y Y_CN b_CN mean_CN

% perform the PLSC and calculate the BSR map
mypls_input;
res  = myPLS_analysis(input,pls_opts);

% the output will directly show how many latent component (LC) is significant
% brain score or latent variable (LV) stores in res.Lx, each coloum is the LV for a LC

BSR_V=res.V./res.boot_results.Vb_std;

% save as gifti file for visiualization
% value of medial wall vertex is 0
BSR_V_full = [MaskData_lh; MaskData_rh]; 
BSR_V_full_LV1(BSR_V_full==1) = bootstrap_ratios_reho(:,1); 
BSR_V_full_LV2(BSR_V_full==1) = bootstrap_ratios_reho(:,2);
BSR_V_full_LV1(BSR_V_full==1,1) = bootstrap_ratios_reho(:,1);
BSR_V_full_LV2(BSR_V_full==1,1) = bootstrap_ratios_reho(:,2);
y_Write(BSR_V_LV1(1:2562, 1),MaskData_lh_gheader,'./lh.LV1.fsaverage4.func.gii');
y_Write(BSR_V_LV1(2563:5124, 1),MaskData_rh_gheader,'./rh.LV1.fsaverage4.func.gii');
y_Write(BSR_V_LV2(1:2562, 1),MaskData_lh_gheader,'./lh.LV2.fsaverage4.func.gii');
y_Write(BSR_V_LV2(2563:5124, 1),MaskData_rh_gheader,'./rh.LV2.fsaverage4.func.gii');

% downsmaple the BSR_V map for decoding

```
