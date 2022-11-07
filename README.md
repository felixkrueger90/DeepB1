# DeepB1

This repository contains Python code for implementing a neural network to estimate complex, channel-wise, relative 2D B1+-maps from a single gradient echo localizer to overcome long subject-specific pTx calibration times as described in [1]. MATLAB algorithms to process the data sets for network traing and subsequent pulse desing are also included. The data sets containing relative B1+ and localizer data of the human body at 7T and the neural network`s weights are available at:https://doi.org/10.6084/m9.figshare.21268650.v1, https://doi.org/10.6084/m9.figshare.21268659.v1, and https://doi.org/10.6084/m9.figshare.21269115.v1. The B1+ maps were computed as described in [2].

##### Authors:
- Felix Krüger                  (<felix.krueger@ptb.de>)
- Christoph S. Aigner           (<christoph.aigner@ptb.de>)
- Kerstin Hammernik             (<k.hammernik@tum.de>)
- Sebastian Dietrich            (<sebastian.dietrich@ptb.de>)
- Max Lutz                      (<max.lutz@ptb.de>)
- Jeanette Schulz-Menger        (<jeanette-schulz-menger@charite.de>)
- Tobias Schaeffter              (<tobias.schaeffter@ptb.de>)
- Sebastian Schmitter           (<sebastian.schmitter@ptb.de>)

Usage
--------

Run the script DataProcessing.m: The user can load the exemplary data set GRE_B1R_109.mat. The data set can be processed to be suitable for network training. Save the variables lvLovalizerInput and lvB1POutput to be used for networkt training. The neural network can be implemented and tranined using NetworkTrainingAndDataPrediction.ipynb. The predicted results can be used for phase-only static pTx pulse design using Pulsedesign.m. The complete network training relies on all data available at https://doi.org/10.6084/m9.figshare.21268650.v1, https://doi.org/10.6084/m9.figshare.21268659.v1, and https://doi.org/10.6084/m9.figshare.21269115.v1.

Contents
--------

##### Test scripts (run these):
    DataProcessing.m                                test script to process the training data 
    NetworkTrainingAndDataPrediction.ipynb          test script to build and train the neural network
    PulseDesign.m                                   test script to evaluate the prediction using phase-only static pTx
##### Routines called by the test scripts:
    b1_phase_shimmingTXfct.m            function to compute subject-specific, static pTx pulses
    catstruct.m                         function to concatenate or merge structures with different fieldnames
    setValues.m                         function to set the parameters for shimming
    quantify_phase_shim_TXfct.m         quantifies the coefficient of variation and the efficiency for phase shimming
    makeColVec.m                        function makes a column vector of the input vector
    Matlab_WS_Startphases_cxXo_1000x64.mat Matrix containing random start phases for shimming 
    
##### Data files used by the test scripts:
    GRE_B1R_109.mat                                               Examplary data set 
    maskHeart.mat                                                 Heart ROI for the examplary data set 
    maskThorax.mat                                                Thorax ROI for the examplary data set 
    
##### External data files used for network traning:
    ProcessedData_loc_96x96x65_b1r_96x96x16_crossValidation.mat          44 data sets for network traning containg channel-wise B1+ maps and complex localizers
    ProcessedData_loc_96x96x65_b1r_96x96x16_testData.mat                 3 data sets for in vivo application containg channel-wise B1+ maps and complex localizers
    
References
------------
[1] Krueger F, Aigner C, Hammernik K, Dietrich S, Lutz M, Menger J, Schaeffter T, Schmitter S. Rapid estimation of 2D relative B1+-maps from localizers in the human heart at 7T using deep learning. Magn. Reson. Med. https://doi.org/10.1002/mrm.29510.

[2] Van de Moortele PF, Ugurbil K. Very Fast Multi Channel B1 Calibration at High Field in the Small Flip Angle Regime. In: Proceedings of the 17th Annual Meeting of ISMRM, Honolulu, USA: Abstract 0367; 2009.

Dependencies
------------
These routines were tested under MATLAB R2020a under Windows, but should also run under older versions.

The data libraries for traning and in vivo application  containing channel-wise B1+ maps and complex localizers, as well as the network`s weights are available at https://doi.org/10.6084/m9.figshare.21268650.v1, https://doi.org/10.6084/m9.figshare.21268659.v1, and https://doi.org/10.6084/m9.figshare.21269115.v1.

The optimization of the static ptX pulse relies on code from Sebastian Schmitter. The function catstruct.m can be found under https://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct, MATLAB Central File Exchange coded by Jos van der Geest.

Please cite appropriately.

