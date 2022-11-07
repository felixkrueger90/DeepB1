%% load data

% load data set containing localizer and corresponding B1+ data for one example
% data set, the whole preprocessed data library for training, cross-validation and in vivo
% application can be found under https://doi.org/10.6084/m9.figshare.21268659 

load('GRE_B1R_109.mat')
load('maskThorax.mat')
load('maskHeart.mat')

%% do shimming

%load start phases
load('Matlab_WS_Startphases_cxX0_1000x64.mat');

% set shim options
optsdef = setValues(cxX0);
lvCounter = 0;

lv_whichSlice = 2;

ROI.masks = permute(maskHeart(:,:,lv_whichSlice),[1 2]);

% calculate shim vector on GT
b1r_chwi_allDataGT = reshape(b1r_chwi_109(:,:,lv_whichSlice,:),[96,96,1,8]).*permute(maskThorax(:,:,lv_whichSlice,:),[1,2,3]);

[betaAll,fAll_GT]     = b1_phase_shimming_TXfct(b1r_chwi_allDataGT*1,ROI.masks,optsdef);
[mmin, mind]    = min(fAll_GT); 
lSolution_GT       = mind;

lSolutionGT = betaAll;
lvShimSetGT = betaAll(:,lSolution_GT);
lvValues = quantify_phase_shim_TXfct(lvShimSetGT,b1r_chwi_allDataGT,ROI.masks); 

% lvShimSetGT is the resulting shim vector, b1r_chwi_allDataGT the b1+ data and lvValues contains the characteristics of the resulting shim 
