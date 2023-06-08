This is the analysis code for "Transdiagnostic neural variability deficits in three major psychiatric disorders"
The preprocessing of the T1w/rest-fMRI is not shown here, which can be easily done with publicly avilable tool like fMRIprep/nilearn/DPABI
Main contents include calculating SDBOLD, performing 4-group partial least squares correlation (PLSC) (https://github.com/MIPLabCH/myPLS) 
then calculating the Bootstraping ratio (BSR) map, and decode the map with neuromaps (https://github.com/netneurolab/neuromaps)

## Step1: calculating the SDBOLD for each individual
This is depending on MATLAB with spm and freesurfer installed

```
%% initializing
data_dir = 'fMRI_surf/'; % data dir contains all the lh* and rh* fmri files
lh_data_list = dir([data_dir 'fsaverage4/lh*']);
rh_data_list = dir([data_dir 'fsaverage4/rh*']);
cd([data_dir 'fsaverage4']);
surf_lh=['lh.pial.fsaverage4.surf.gii'];
surf_rh=['rh.pial.fsaverage4.surf.gii'];

% get mask
MaskData_lh_gheader=gifti('lh.cortex.fsaverage4.func.gii');
MaskData_rh_gheader=gifti('rh.cortex.fsaverage4.func.gii');
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

##step 3 decoding the disruption pattern
the gene-PLSC has a fully documented pipeline in https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md
the symptom-PLSC can also be conducted 
so here we only show the pipeline of decoding the pattern with neurosynth cognitive terms and neurotransmiiter map (neuromaps)

notice that this pipeline is performed under python

```
# import needed lib
import contextlib, json, os, warnings, requests, nibabel, abagen
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
import neuromaps
from neuromaps.datasets import available_annotations
import numpy as np
from netneurotools import datasets as nntdata
from neuromaps.parcellate import Parcellater
from neuromaps.images import annot_to_gifti
from neuromaps import datasets
from neuromaps import transforms
from neuromaps.parcellate import Parcellater
from neuromaps.images import annot_to_gifti
from neuromaps import images, nulls
from neuromaps import stats
import scipy.io as sio
from pathlib import Path
import pandas as pd
import sys
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

schaefer = nntdata.fetch_schaefer2018('fsaverage5')['400Parcels7Networks']
parc_gifti = annot_to_gifti(schaefer)
parc_fsaverage5 = Parcellater(parc_gifti, 'fsaverage')
parc_mni = Parcellater('/home/weiw/Desktop/data/MicroStruProfile/aparc.2009s.MNIvolumeSpace/test_HCP_MMP1.nii.gz', 'mni152')

expression = abagen.get_expression_data(parc_gifti, ibf_threshold=0.5, tolerance=2, sim_threshold=2, missing='interpolate', n_proc=6)
expression.to_csv('/home/weiw/weiw_data_hdd/ABA_Schaefer_7_400_data_fsaverage.csv')

Ind=list(expression.index)
Col=list(expression.columns)
sio.savemat('/home/weiw/weiw_data_hdd/neuromaps_data/predictor/abagen_expression_400schaefer17_rnaseq.mat', {'data_expression':expression.values, 'index':Ind, 'cols_expression':Col})

data_PET = datasets.fetch_annotation(tags='PET')
data_meta = datasets.fetch_annotation(source='raichle')

parc_data = np.loadtxt('/home/weiw/weiw_data_hdd/MicroStruProfile/BSR_gradient_1_LV1.txt')
for keys in data_PET:
  if 'beliveau2017' in keys:
    if 'MNI152' in keys:
      continue
    else:
      fsaverage5_data = transforms.fsaverage_to_fsaverage(data_PET[keys], '10k')
      test_parc = parc_fsaverage5.fit_transform(fsaverage5_data, 'fsaverage')
      parc_data = np.column_stack((parc_data, test_parc.T))
  elif 'norgaard2021' in keys:
    if 'MNI152' in keys:
      continue
    else:
      fsaverage5_data = transforms.fsaverage_to_fsaverage(data_PET[keys], '10k')
      test_parc = parc_fsaverage5.fit_transform(fsaverage5_data, 'fsaverage')
      parc_data = np.column_stack((parc_data, test_parc.T))
  else:
    test_parc = parc_mni.fit_transform(data_PET[keys], 'mni152')
    parc_data = np.column_stack((parc_data, test_parc.T))

for keys in data_meta:
  fsaverage5_data = transforms.fslr_to_fsaverage(data_meta[keys], '10k')
  test_parc = parc_fsaverage5.fit_transform(fsaverage5_data, 'fsaverage')
  parc_data = np.column_stack((parc_data, test_parc.T))

# get roatated table
rotated1 = nulls.burt2020(parc_data[:,0], atlas='fsaverage', density='10k', n_perm=1000, seed=7, parcellation=parc_gifti, n_proc=20)

#rotated1 = nulls.alexander_bloch(parc_data[:,1], atlas='fsaverage', density='10k', n_perm=10000, seed=7, parcellation=parc_gifti)

for i in range(3, 44):
 corr, pval = stats.compare_images(parc_data[:,0], parc_data[:,i], nulls=rotated1)
 print(f'{i-2} r = {corr:.3f}, p = {pval:.4f}')
 
rotated2 = nulls.alexander_bloch(parc_data_norm[:,1], atlas='fsaverage', density='10k', n_perm=1000, seed=7, parcellation=parc_gifti, n_proc=20)
for i in range(3, 44):
 corr, pval = stats.compare_images(parc_data[:,1], parc_data[:,i], nulls=rotated2)
 print(f'{i-2} r = {corr:.3f}, p = {pval:.4f}')
 
rotated3 = nulls.burt2020(parc_data[:,2], atlas='fsaverage', density='10k', n_perm=1000, seed=7, parcellation=parc_gifti, n_proc=20)
for i in range(0, 44):
 corr, pval = stats.compare_images(parc_data[:,2], parc_data[:,i], nulls=rotated3)
 print(f'{i-2} r = {corr:.3f}, p = {pval:.4f}')
 

parc_data = np.loadtxt('/home/weiw/weiw_data_hdd/MicroStruProfile/T_gradient_1_reg.txt')
#parc_data = preprocessing.scale(parc_data)
rotated1 = nulls.burt2020(parc_data[:,0], atlas='fsaverage', density='10k', n_perm=1000, seed=777, parcellation=parc_gifti, n_proc=20)
rotated2 = nulls.burt2020(parc_data[:,1], atlas='fsaverage', density='10k', n_perm=1000, seed=777, parcellation=parc_gifti, n_proc=20)
rotated3 = nulls.burt2020(parc_data[:,2], atlas='fsaverage', density='10k', n_perm=1000, seed=777, parcellation=parc_gifti, n_proc=20)

sio.savemat('/home/weiw/Desktop/data/MicroStruProfile/young_nullmaps.mat', {'young_nullmaps':rotated1})
sio.savemat('/home/weiw/Desktop/data/MicroStruProfile/mean_nullmaps.mat', {'mean_nullmaps':rotated2})
sio.savemat('/home/weiw/Desktop/data/MicroStruProfile/old_nullmaps.mat', {'old_nullmaps':rotated3})

# get cognitive atlas
url = 'https://cognitiveatlas.org/api/v-alpha/concept'
req = requests.get(url)
req.raise_for_status()
concepts = set([f.get('name') for f in json.loads(req.content)])

# get neurosynth terms
ns_terms=pd.read_csv('/home/weiw/weiw_data_hdd/neuromaps_data/predictor/neurosynth_terms.csv')
ns_terms=np.array(ns_terms)
ns_terms_set=set()
for item in ns_terms:
 for i in item:
  ns_terms_set.add(i)

# load&parcellate neurosynth meta map
parc_ns_data = np.loadtxt('/home/weiw/weiw_data_hdd/MicroStruProfile/T_gradient_1_reg.txt')

usable_terms=ns_terms_set&concepts
ns_dataDIR='/home/weiw/neurosynth_data/raw/brainstat_neurosynth/'
usable_terms=list(usable_terms)
for i in range(0, len(usable_terms)):
 ns_map=images.load_nifti(ns_dataDIR + 'Neurosynth_TFIDF__' + usable_terms[i] + '_z_desc-consistency.nii.gz')
 parc_temp= parc_mni.fit_transform(ns_map, 'mni152')
 parc_ns_data = np.column_stack((parc_ns_data, parc_temp.T))
 

parc_ns_data_norm = scale_values(parc_ns_data, 0, 100, axis=0)
#np.savetxt('/home/weiw/Desktop/data/neuromaps_data/predictor/neurosynth_400schaefer17neworks.csv', parc_ns_data, delimiter=',')

orig_stdout = sys.stdout
f = open('young_neurosynth.csv', 'w')
sys.stdout = f
for i in range(3, len(usable_terms)+3):
 corr, pval = stats.compare_images(parc_ns_data[:,0], parc_ns_data[:,i], nulls=rotated1)
 print(f'{i-2} {usable_terms[i-3]}, r = {corr:.3f}, p = {pval:.4f}')
 
sys.stdout= orig_stdout
f.close() 
 
orig_stdout = sys.stdout
f = open('old_neurosynth.csv', 'w')
sys.stdout = f
for i in range(3, len(usable_terms)+3):
 corr, pval = stats.compare_images(parc_ns_data[:,2], parc_ns_data[:,i], nulls=rotated3)
 print(f'{i-2} {usable_terms[i-3]}, r = {corr:.3f}, p = {pval:.4f}')

sys.stdout= orig_stdout
f.close() 




```
