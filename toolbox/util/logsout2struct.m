function out = logsout2struct(LOGSOUT,OutputType)
%% Convert logsout to structure

% Input processing
narginchk(1,2);
nin = nargin;
if nin < 2
    OutputType = 'double';
end

% Return output of the requested type in a structure
sigNames =  getElementNames(LOGSOUT);
out = struct([]);
switch OutputType    
    case 'double'
        for i = 1:length(sigNames)
            tempSignal = get(LOGSOUT,sigNames{i});
            if isequal(i,1)
                out(1).Time = tempSignal.Values.Time;                
            end
            out(1).(sigNames{i}) = tempSignal.Values.Data;
        end
        
    case 'tvmat'        
        for i = 1:length(sigNames)
            eval([sigNames{i} '= signal2tvmat(get(LOGSOUT,sigNames{i}));']);
            out(1).(sigNames{i}) = eval(sigNames{i});
        end
        
    case 'ts'
        for i = 1:length(sigNames)
            tempSignal = get(LOGSOUT,sigNames{i});
            out(1).(sigNames{i}) = tempSignal.Values;
        end
        
    otherwise
        error('Invalid Option');
end