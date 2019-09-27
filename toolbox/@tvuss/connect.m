function sysc = connect(varargin)

%% Split apart systems
% Syntax from connect documentation:
%  1) sysc = connect(sys1,...,sysN,inputs,outputs)
%  2) sysc = connect(sys1,...,sysN,inputs,outputs,APs)
%  3) sysc = connect(blksys,connections,inputs,outputs)
%  4) sysc = connect(___,opts)
% #1 and #2 have systems followed by character strings to specify I/O/APs.
% #3 has systems followed by double arrays to specify I/Os.
% Options are specified by either Name/Value pairs (with Name as a char)
% or a connectOptions object.
nin = nargin;
N = nin;
for i=1:nin
    vi = varargin{i};
    if ischar(vi) || iscellstr(vi) || isa(vi,'ltioptions.connect')
        N = i-1;
        break;
    end
end

%% Use switchyard to handle object interactions
[isTimeInvariant,varargin{1:N}] = tvswitchyard( varargin{1:N} );

%% Handle different TVMAT data types
if isTimeInvariant
    for k=1:N
        varargin{k} = varargin{k}.Data;
    end
    Data = connect( varargin{:} );
    
    if isa(Data,'ss') || isa(Data,'tf') || isa(Data,'zpk')
        sysc = tvss(Data);
    elseif isa(Data,'uss')
        sysc = tvuss(Data);
    end
    
else
    % Connect works on system arrays 
    Time = varargin{1}.Time;
    IM = varargin{1}.InterpolationMethod;
    for k=1:N
        varargin{k} = varargin{k}.Data;
    end
    CData = connect( varargin{:} );
    
    if isa(CData,'ss') || isa(CData,'tf') || isa(CData,'zpk')
        sysc = tvss(CData,Time,IM);
    elseif isa(CData,'uss')
        sysc = tvuss(CData,Time,IM);
    end
end


