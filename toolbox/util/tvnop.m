function varargout = tvnop(opfh,varargin)
% TVNOP N-ary operations for a TVMAT
%
% TVNOP is a utility function that handles N-ary operations (e.g. eig, svd,
% etc.) on TVMATs. The user should not call TVNOP directly. Instead each
% N-ary operation has a front end function that the user calls. The
% front end function then call calls TVNOP to do the required calculations.

%% Input / output error checking
nout = nargout;
varargout = cell(1,nout);

%% Separate char inputs from remainder of data
idx = cellfun(@ischar,varargin);
optstr = varargin(idx);
varargin = varargin(~idx);
nin = numel(varargin);

%% Use switchyard to handle object interactions
[isTimeInvariant,varargin{:}] = tvswitchyard(varargin{:});

%% Handle different TVMAT data types
if isTimeInvariant
    for k=1:nin
        varargin{k} = varargin{k}.Data;
    end
    [varargout{1:nout}] = opfh(varargin{:},optstr{:});
    
    for k=1:nout
        CData = varargout{k};
        if isa(CData,'double') || isa(CData,'logical')
            varargout{k} = tvmat(CData);
        elseif isa(CData,'ss') || isa(CData,'tf') || isa(CData,'zpk')
            varargout{k} = tvss(CData);
        elseif isa(CData,'umat')
            varargout{k} = tvumat(CData);
        elseif isa(CData,'uss')
            varargout{k} = tvuss(CData);
        end
        
        %        varargout{k} = tvmat( varargout{k} );
    end
else
    Time = varargin{1}.Time;
    IM = varargin{1}.InterpolationMethod;
    Ts = varargin{1}.Ts;
    Nt = numel(Time);
    vinI = cell(1,nin);
    voutI = cell(1,nout);
    voutALL = cell(1,nout);
    for i=1:Nt
        % Grab Data associated with t(i) from all inputs
        for k=1:nin
            vinI{k} = varargin{k}.Data(:,:,i);
        end
        
        % Perfom Operation
        [voutI{:}] = opfh(vinI{:},optstr{:});
        
        % Store Data associated with operation performed at t(i)
        if i==1
            for k=1:nout
                voutALL{k} = zeros([numel(voutI{k}) Nt]);
            end
        end
        for k=1:nout
            voutALL{k}(:,i) = voutI{k}(:);
        end
    end
    
    % Reshape outputs
    if isequal(Ts,0)
        Arg3 = IM;        
    else
        Arg3 = Ts;
    end    
    for k=1:nout
        CData = reshape(voutALL{k},[size(voutI{k}) Nt]);
        if isa(CData,'double') || isa(CData,'logical')                            
            varargout{k} = tvmat(CData,Time,Arg3);            
        elseif isa(CData,'ss') || isa(CData,'tf') || isa(CData,'zpk')
            varargout{k} = tvss(CData,Time,IM);
        elseif isa(CData,'umat')
            varargout{k} = tvumat(CData,Time,IM);
        elseif isa(CData,'uss')
            varargout{k} = tvuss(CData,Time,IM);
        end
        
        %        voutALL{k} = reshape(voutALL{k},[size(voutI{k}) Nt]);
        %        varargout{k} = tvmat(voutALL{k},Time,IM);
    end
end

