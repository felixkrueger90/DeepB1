%% load data

% this is an example for processing the data for the neural network, the 
% whole preprocessed data set can be found under  https://doi.org/10.6084/m9.figshare.21268659 

%load input localizer and  groundtruth B1P for one example data set
load('GRE_B1R_109.mat')
% load thorax ROI
load('maskThorax.mat')

b1r_chwi_allData = b1r_chwi_109;
gre_chwi_allData = gre_chwi_109;
clear b1r_chwi_109 gre_chwi_109

szSlectB1Rp = size(b1r_chwi_allData,3);
szNoChannels = size(b1r_chwi_allData,4);

%% prepare data for network

lvCounter = 0;

for m = 1:szSlectB1Rp
    clear b1r_chwi_CurrSlice gre_chwi_CurrSlice currentMask xCurrentMask yCurrentMask gre_chwi_CurrSliceAll  b1r_chwi_CurrSliceAll;
    b1r_chwi_CurrSlice1 = b1r_chwi_allData(:,:,m,:);
    b1r_chwi_CurrSlice1 = (reshape(b1r_chwi_CurrSlice1,[96,96,8]));
     
    szB1rchwiCurrSlice = size(b1r_chwi_CurrSlice1);
    
    gre_chwi_CurrSlice1 = gre_chwi_allData(:,:,m,:);
    gre_chwi_CurrSlice1 = (reshape(gre_chwi_CurrSlice1,[96,96,32]));
       
    szGre_chwi_CurrSlice = size(gre_chwi_CurrSlice1);
    
    clear currMaskThorax;
    currMaskThorax = maskThorax(:,:,m);

    b1r_chwi_CurrSlice = b1r_chwi_CurrSlice1;
    
    gre_chwi_CurrSlice = gre_chwi_CurrSlice1;
    
    b1r_chwi_CurrSliceReal = (real(b1r_chwi_CurrSlice)).*currMaskThorax;
    b1r_chwi_CurrSliceImag = (imag(b1r_chwi_CurrSlice)).*currMaskThorax;
    
    b1r_chwi_CurrSliceAll(:,:,1:2:16) = b1r_chwi_CurrSliceReal;
    b1r_chwi_CurrSliceAll(:,:,2:2:16) = b1r_chwi_CurrSliceImag;
    
    slctChWiseGREAugmSum = gre_chwi_CurrSlice;
    b1rSquare = abs(slctChWiseGREAugmSum);
    lvLovalizerAbs = rescale(sqrt(sum(b1rSquare,3)));
    
    RxReal = (real(slctChWiseGREAugmSum)).*currMaskThorax;
    RxImag = (imag(slctChWiseGREAugmSum)).*currMaskThorax;
    
    gre_chwi_CurrSliceAll(:,:,1) = lvLovalizerAbs.*currMaskThorax;
    gre_chwi_CurrSliceAll(:,:,2:2:65) = RxReal;
    gre_chwi_CurrSliceAll(:,:,3:2:65) = RxImag;

    slctChWiseB1R = b1r_chwi_CurrSliceAll;
                
    slctChWiseGRE = gre_chwi_CurrSliceAll;


    lvCounter = lvCounter+1;

    b1r_chwi_CurrSliceMax = max(max(max(abs(slctChWiseB1R(:,:,1:2:16)+1i*slctChWiseB1R(:,:,2:2:16)))));

    lvB1POutput(:,:,:,lvCounter) = single(slctChWiseB1R./b1r_chwi_CurrSliceMax); %norm on max value

    gre_chwi_CurrSliceMax = max(max(max(abs(slctChWiseGRE(:,:,2:2:65)+1i*slctChWiseGRE(:,:,3:2:65)))));

    lvLovalizerInput(:,:,1,lvCounter) = single(slctChWiseGRE(:,:,1)); %norm on max value

    lvLovalizerInput(:,:,2:65,lvCounter) = single(slctChWiseGRE(:,:,2:65)./gre_chwi_CurrSliceMax);
  
end 

% lvLovalizerInput is utilized as the input and lvB1POutput as the output
% of the neural network

