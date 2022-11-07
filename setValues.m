function optsdef = setValues(cxX0)

opts_B1shim.ALGO = 'homogeneity';
% opts_B1shim.ALGO = 'efficiency';

% opts_B1shim.ALGO = 'forceefficiencytovalue';

optsdef.BETA0 = cxX0;
optsdef.ALGO = opts_B1shim.ALGO;
optsdef.MASKSTAYSAME = [];
optsdef.SUMUPTOTHATTHRESHOLD = 0.1;
optsdef.NOOFSTARTPHASES = 32;
optsdef.EFFVALUEFORCE = 0.5;
optsdef.MEANPOWER = 1; 
optsdef.LAMBDA = 0;
