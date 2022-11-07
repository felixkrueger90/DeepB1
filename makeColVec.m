function colVec = makeColVec(inVec)
%% ************************************************************************
%
%   makes a column vector of the input vector inVec: does not matter if it
%   is alreay a col vec or not. Does  work with matrices, too, but in this
%   case should be handled with care!
%
%
%   INPUT:                                                  [unit]
%   ----------------------------------------------------------------
%   inVec:      input vector
%
%   OUTPUT:
%   ----------------------------------------------------------------
%   colVec  	output column vector
%
%%************************************************************************


    vecsize = size(inVec);
    
    if(max(size(vecsize)) > 2 || (min(vecsize) > 1))
        colVec = inVec;
        fprintf('makeColVec: use vectors only!\n');
        return ;
    end
    
    if(vecsize(2) > vecsize(1))
        colVec = inVec.';
    else
        colVec = inVec;
    end

end