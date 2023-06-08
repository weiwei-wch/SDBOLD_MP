This is the analysis code for "Transdiagnostic neural variability deficits in three major psychiatric disorders"
The preprocessing of the T1w/rest-fMRI is not shown here, which can be easily done with publicly avilable tool like fMRIprep/nilearn/DPABI
Main contents include calculating SDBOLD, performing 4-group partial least squares correlation (PLSC) (https://github.com/MIPLabCH/myPLS) 
then calculating the Bootstraping ratio (BSR) map, and decode the map with neuromaps (https://github.com/netneurolab/neuromaps)

## calculating the SDBOLD for each individual
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
MaskData_lh_gheader=gifti('/home/weiwei/Desktop/polar_LDA/lh.cortex.fsaverage4.func.gii');
MaskData_rh_gheader=gifti('/home/weiwei/Desktop/polar_LDA/rh.cortex.fsaverage4.func.gii');
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

sd_data = [sd_lh_z, sd_rh_z]
```
