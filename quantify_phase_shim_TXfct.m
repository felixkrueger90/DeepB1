function  ValueStruct = quantify_phase_shim_TXfct(shimset, b1p, roi, varargin)
 %% ************************************************************************
%
%   Quantifies the coefficient of variation and the efficiency for phase
%   shimming
%   
%   Author: S.Schmitter
%   Date: Jan 2016
%
%
%   dependencies:
%   - multiprod
%   - parseVariableInputs
%   - catstruct
%
%   INPUT:                                                  [unit]
%   ----------------------------------------------------------------
%   shimset     coplex vector of phase shim values 
%   b1p         b1+ maps
%   roi         region of interest
%
%   Options (with standard prefs)                           [unit]
%   ----------------------------------------------------------------
%                          
%
%   OUTPUT:
%   ----------------------------------------------------------------
%   ValueStruct     Struct containing the quantified values 
%
%% ************************************************************************

    optsdef.TEST = 0;
    optsdef.SAVEMAPS = 0;

    opts    = catstruct(optsdef,parseVariableInputs(varargin));

    shimvec = makeColVec(shimset);
    
    b1pat = multiprod(b1p,shimvec,4,1);
    b1sumofmag = multiprod(abs(b1p),abs(shimvec),4,1);
    
    EfficiencyMap = abs(b1pat)./b1sumofmag;
    ValueStruct.Efficiency = mean(EfficiencyMap(~~roi));
    ValueStruct.EfficiencyMin = min(EfficiencyMap(~~roi));
    ValueStruct.EfficiencyMax = max(EfficiencyMap(~~roi));
    
    tmp = abs(b1pat(~~roi));
    
    ValueStruct.CV = std(tmp)/mean(tmp);
    
    if(opts.SAVEMAPS)
        ValueStruct.B1map = b1pat;
        ValueStruct.EfficiencyMap = EfficiencyMap;
    else
        ValueStruct.B1map = [];
        ValueStruct.EfficiencyMap = [];
    end
    ValueStruct.ROIsize = sum(roi(:));
end