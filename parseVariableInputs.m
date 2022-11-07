function prop = parseVariableInputs(arguments)
%% **********************************************************************
%
% function prop = parseVariableInputs(varargin)
%
%   usage: prop = parseVariableInputs({'MYPARAM1 = 2';'MYPARAM2 =
%   4'},'mySpecialFilenameWithoutExtension','My general comment');
%
%   Autor:  Sebastian Schmitter
%   Date:   Sept. 2011
% 
%
%   INPUT:                                                          [unit]
%   -----------------------------------------------------------------------
%   first parameter:
%   the varargin-parameters from any functions in cell format
%   e.g. {  'BWTP = 8';...
%           'ISHANNING = true';...
%           'DURATION = 1024E-6';...
%           'THICKNESS = 5E-3'}
%
%   second parameter:
%   filename for saving etc.; will be stored in prop.FILENAME
%
%   third parameter:
%   arbitrary comment; will be stored in prop.COMMENT
%
%   OUTPUT:
%   -----------------------------------------------------------------------
%   properties struct, e.g.
%
%   prop.BWTP 
%   prop.ISHANNING
%   etc.
%           
%% **********************************************************************

    optargin = size(arguments,2);

    %first check wheater prop is a struct or a cell
    if(optargin > 0)
        if(iscell(arguments))
           z = arguments{1,1};
           if(isstruct(z))
               prop = z;
           else
               for lL=1:length(arguments{1,1})
                   eval(['prop.',char(arguments{1,1}(lL)),';']);
               end
           end
        else %if iscell
        %otherwise it is assumed that inprop is a struct
            error('varargin must have cell format');
        end
        if(optargin >1)
            %second input parameter
            prop.FILENAME   = char(arguments{1,2});
            
        end
        if(optargin >2)
            prop.COMMENT    =  char(arguments{1,3});
        end
        
    else
        prop    = struct;
    end %if optargin

end