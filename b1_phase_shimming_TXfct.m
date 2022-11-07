function [betaAll,fAll] = b1_phase_shimming_TXfct(b1p,roi,varargin)

%% ************************************************************************
%
% This is the optimization routine for predistorted B1shim. It uses a
% constraint nonliniear optimization algorithm to vary the phases (only) in
% order to match b1 plus with the target. 
% 
%
%   INPUT:                                                  [unit]
%   ----------------------------------------------------------------
%   b1plusmph:  complex 4D data with b1 plus profiles 
%   dMask:      contains the mask for the optimization ROI. The dimension 
%               of the mask must be 3d with same size as dTarget as well as
%               the size of the first 3 dim of b1plusmph
%   dTarget:    contains the target in a.u.
%
%
%   OUTPUT:
%   ----------------------------------------------------------------
%   xhat:       complex vector with phase setting of the cannels
%
%% ************************************************************************
% opts    = catstruct(parseVariableInputs(varargin));
% optsdef.ALGO   = 'HOMOGENEITY';
% optsdef.LAMBDA = 0; %i.e. pure  homogeneity constraint
optsdef.WHICHCHANNELS = (1:size(b1p,4)).';
% optsdef.MASKSTAYSAME = [];
% optsdef.BETA0 = ones(size(b1p,4),1);  %this is the starting b1 shim value
% optsdef.NOOFSTARTPHASES = 1;
% optsdef.MEANPOWER = 1;
% optsdef.SUMUPTOTHATTHRESHOLD = 0.1;  %what is that?
% optsdef.EFFVALUEFORCE = 0.5;


optargin = size(varargin,2);

opts    = catstruct(optsdef,parseVariableInputs(varargin));

opts.CHANNELVEC = [opts.WHICHCHANNELS; opts.WHICHCHANNELS];%+size(b1plusmph,4)];

%shim vector (starting value
Beta0 = opts.BETA0;

b1p	= double(b1p);        %not clear if that is needed..
b1p(~isfinite(b1p))     = 0;    %erases inf. values

scfact      = 1;%/mean(abs(b1p(:)));%scaling factor, needed later

b1psize     = size(b1p);
sz = b1psize;
lNoOfCha	= b1psize(4);

b1plustmp	= reshape(b1p,[prod(sz(1:3)),sz(4)]);
Atmp        = b1plustmp(find(roi),:); %optimize only those values withinn mask
                                            %creates Nx16 matrix (for 16
                                            %channels)
Ttmp        = roi(find(roi));     %shrink the target to 1D array
%find nonzero values, i.e. those for A*ones(16,1) != 0
lNonZeroInd = find(Atmp*ones(lNoOfCha,1));
A           = scfact*Atmp(lNonZeroInd,:);
T           = Ttmp(lNonZeroInd);

% mod ss: if the overall pattern should not change:
B           = b1p(find(opts.MASKSTAYSAME),:);
opt.B      = B;

opts.LASTLOWINDEX = ceil(opts.SUMUPTOTHATTHRESHOLD.*size(A,1));

%initial point setting:
%cxX0        = ones(lNoOfCha,1)/2;
%in case of randomized initial point setting, use:
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)))
%cxX0 = exp(2*pi*1i*rand(lNoOfCha,1,'double'));

betaAll = zeros(sz(4), opts.NOOFSTARTPHASES);
fAll    = zeros(opts.NOOFSTARTPHASES,1);

%run loop over different starting values:
% parfor lL = 1:opts.NOOFSTARTPHASES
    for lL = 1:opts.NOOFSTARTPHASES
    
    cxBeta0 = makeColVec(opts.BETA0(lL,1:sz(4)));
    dBeta0         = [real(cxBeta0); imag(cxBeta0)]; %only real values, but 2*lNoOfCha


    %call of kernel (fmincon routine):
    [x,fval,exitflag,output] = runOptimization(A,dBeta0,T,opts);

    betaAll(:,lL)       = makeColVec(x(1:end/2)+1i*x(end/2+1:end)); 
    fAll(lL,1)      	= fval;
end
%end run loop


end



%% *************************************************************************
% local Functions:
%%*************************************************************************

function [x,fval,exitflag,output] = runOptimization(A,dX0,dTVec,prop)
    %output initial coefficient value:
    
    
    %options:
%     options = optimset('Algorithm','interior-point');
%     options = optimset('Algorithm','interior-point','display','none');
    options = optimoptions(@fmincon,'Algorithm','interior-point');
% options = optimoptions(@fmincon,'Algorithm','interior-point');
    %options.LargeScale      = 'off';
    options.MaxIter         = 3e3;
%     options.LargeScale      = 'on';
    options.MaxFunEvals     = 1e5;
    options.TolX            = 1e-12;
    options.TolCon          = 1e-12;
    options.TolFun          = 1e-12;
    
    meanpower = prop.MEANPOWER;
    
    lChaVec = prop.CHANNELVEC;
    
    switch lower(prop.ALGO)
        case 'homogeneity' 
            [x,fval,exitflag,output]  = fmincon(@(x)objfunHomogeneity(x,A,dTVec,meanpower,prop),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
       	case 'homogeneityold' 
            [x,fval,exitflag,output]  = fmincon(@(x)objfunHomogeneityOld(x,A,dTVec,meanpower),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
        case 'efficiency'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunEfficiency(x,A,dTVec,meanpower,prop),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
        case 'efficiencyold'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunEfficiencyOld(x,A,dTVec,meanpower),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
        case 'effhomhybrid'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunEffHomHybrid(x,A,dTVec,meanpower,prop.LAMBDA),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);      
        case 'effhomhybridtikho'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunEffHomHybridTikho(x,A,dTVec,meanpower,prop.LAMBDA),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);     
                                    
       case 'homenergytikho'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunHomEnergyTikho(x,A,dTVec,meanpower,prop.LAMBDA),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);    
                                    
        case 'homogenefficiency'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunHomogenEfficiency(x,A,dTVec),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
        case 'forceefficiencytovalue'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunForceEfficiencyToValue(x,A,dTVec,prop.EFFVALUEFORCE),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
        case 'shimblackholeonly'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunBlackHoleShim1(x,A,dTVec,meanpower,prop),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
        case 'shimblackhole2'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunBlackHoleShim2(x,A,dTVec,meanpower,prop),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);
                                                                
        case 'oneovermin'
            [x,fval,exitflag,output]  = fmincon(@(x)objfunOneOverMin(x,A,dTVec,meanpower,prop),...
                                        dX0,...
                                        [],[],[],[],[],[],...
                                        @confun,...
                                        options);                                    
                                    
    end
    
end


    
%objective function:
function f =  objfunOneOverMin(xx,A,dTVec,meanpower,prop)
    
     xtmp = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    %x = makeColVec(prop.BETA0);
    %x(prop.WHICHCHANNELS) = xtmp;
    %atmp        = abs(A*x);
    %sca         = max(size(atmp));
    x = xtmp;
    
    pat = abs(A*x);
    
    f        = -(min(pat(:)));
    
%     meana       = mean(atmp);
%     meant       = mean(dTVec);
%     stda        = std(atmp-dTVec/meant*meana);
%     fhom        = stda/(meana)^meanpower;   
%     
    %f           = feff;
 
end
    
%objective function:
function f = objfunBlackHoleShim2(xx,A,dTVec,meanpower,prop)
    
    xtmp = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    x = prop.BETA0;
    x(prop.WHICHCHANNELS) = xtmp;
   
    %this is the prediciton for the Black hole area
    atmp        = abs(A*x);
    %this the prediciton for the non-BH area (that should stay constant)
    btmp        = abs(prop.B*x);
    b0          = abs(prop.B*prop.X0);

    %sca         = max(size(atmp));
    
     feff        = -sum(abs(A*x)./(abs(A)*abs(x)));
     sb         = sum((btmp-b0).^2);
%      
%     meana       = mean(atmp);
%     meant       = mean(dTVec);
%     stda        = std(atmp-dTVec/meant*meana);
%     fhom        = stda/(meana)^meanpower;   
    
    f           = (sb+prop.LAMBDA*feff);%+feff/sca/8;
 
end




%objective function:
function f = objfunBlackHoleShim1(xx,A,dTVec,meanpower,prop)
    
    xtmp = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    x = prop.BETA0;
    x(prop.WHICHCHANNELS) = xtmp;
    
    atmp        = abs(A*x);
    
    aa          = sort(atmp);
    amed        = median(aa);
    afin        = sum((amed-aa(1:ceil(prop.LASTLOWINDEX))).^2);
    %sca         = max(size(atmp));
    
%     %feff        = -sum(abs(A*x)./(abs(A)*abs(x)));
%     
%     meana       = mean(atmp);
%     meant       = mean(dTVec);
%     stda        = std(atmp-dTVec/meant*meana);
%     fhom        = stda/(meana)^meanpower;   
    
    f           = (afin);%+feff/sca/8;
  
end


function f = objfunHomogeneity(xx,A,dTVec,meanpower,prop)
    
    xtmp = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    %these two lines use if only selected channels should be optimized.
    %does not correctlywork at the omeent
    x = prop.BETA0;
    %x(prop.WHICHCHANNELS) = xtmp;
    
    
    x = xtmp;
    %x = x.';
    atmp        = abs(A*x);
    %sca         = max(size(atmp));
    
    %feff        = -sum(abs(A*x)./(abs(A)*abs(x)));
    
    meana       = mean(atmp);
    meant       = mean(dTVec);
    stda        = std(atmp-dTVec/meant*meana);
    fhom        = stda/(meana)^meanpower;   
    
%     atmp2        = angle(A*x);
%     %sca         = max(size(atmp));
%     
%     %feff        = -sum(abs(A*x)./(abs(A)*abs(x)));
%     
%     meana2       = mean(atmp2);
%     meant2       = mean(dTVec);
%     stda        = std(atmp2-dTVec/meant2*meana2);
%     fhom2        = stda/(meana2)^meanpower;   
    
    
    f           = fhom;%+fhom2;%+feff/sca/8;
end

function f = objfunHomogeneityOld(xx,A,dTVec,meanpower,prop)
    
    xtmp = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    x = xtmp;%prop.BETA0;
    %x(prop.WHICHCHANNELS) = xtmp;
    
    atmp        = abs(A*x);
    %sca         = max(size(atmp));
    
    %feff        = -sum(abs(A*x)./(abs(A)*abs(x)));
    
    meana       = mean(atmp);
    meant       = mean(dTVec);
    stda        = std(atmp-dTVec/meant*meana);
    fhom        = stda/(meana)^meanpower;   
    
    f           = fhom;%+feff/sca/8;
end

%objective function:
function f = objfunEfficiency(xx,A,dTVec,meanpower,prop)
    xtmp = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    x = xtmp;%x = (prop.BETA0);
    %x(prop.WHICHCHANNELS) = xtmp;
    %atmp        = abs(A*x);
    %sca         = max(size(atmp));
    
    feff        = -sum(abs(A*x)./(abs(A)*abs(x)));
    
%     meana       = mean(atmp);
%     meant       = mean(dTVec);
%     stda        = std(atmp-dTVec/meant*meana);
%     fhom        = stda/(meana)^meanpower;   
%     
    f           = feff;
end

%objective function:
function f = objfunEfficiencyOld(xx,A,dTVec,meanpower,prop)
    x = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    %atmp        = abs(A*x);
    %sca         = max(size(atmp));
    
    feff        = -sum(abs(A*x)./(abs(A)*abs(x)));
    
%     meana       = mean(atmp);
%     meant       = mean(dTVec);
%     stda        = std(atmp-dTVec/meant*meana);
%     fhom        = stda/(meana)^meanpower;   
%     
    f           = feff;
end





%objective function:
function f = objfunForceEfficiencyToValue(xx,A,dTVec,effvalue)
    x = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));


    
    feff        = (abs(A*x)./(abs(A)*abs(x)));
    
    stdeff        = sum((feff(:)-effvalue).^2);

    f           = stdeff;
end

%objective function:
function f = objfunEffHomHybrid(xx,A,dTVec,meanpower,dLambda)
    x = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    
    atmp        = abs(A*x);
    %sca         = max(size(atmp));
    
    feff        = mean(abs(A*x)./(abs(A)*abs(x)));
    
     meana       = mean(atmp);
     meant       = mean(dTVec);
     stda        = std(atmp-dTVec/meant*meana);
     fhom        = stda/(meana)^meanpower;   
%     
    f           = (1-dLambda)*fhom*fhom+dLambda*1/feff/feff;
end

%objective function:
function f = objfunEffHomHybridTikho(xx,A,dTVec,meanpower,dLambda)
    x = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    
    atmp        = abs(A*x);
    %sca         = max(size(atmp));
    
    feff        = mean(abs(A*x)./(abs(A)*abs(x)));
    
     meana       = mean(atmp);
     meant       = mean(dTVec);
     stda        = std(atmp-dTVec/meant*meana);
     fhom        = stda/(meana)^meanpower;   
%     
    f           = fhom*fhom+dLambda^2*1/feff/feff;
end

%objective function:
function f = objfunHomEnergyTikho(xx,A,dTVec,meanpower,dLambda)
    x = makeColVec(xx(1:end/2)+1i*xx(end/2+1:end));

    
    atmp        = abs(A*x);
    %sca         = max(size(atmp));
    
    %feff        = mean(abs(A*x)./(abs(A)*abs(x)));
    
    meana       = mean(atmp(:));
    meant       = mean(dTVec(:));
    stda        = std(atmp/meana*meant-dTVec);
    %fhom        = stda/(meana)^meanpower;   
%     
    f           = stda.^2+dLambda^2*sum(abs(x*meant/meana).^2);
end

% function for constraints 
function [cineq, ceq] = confun(xx)
    x       = xx(1:end/2)+1i*xx(end/2+1:end); 
    ceq     = [abs(x)-1];
    cineq   = [];
end