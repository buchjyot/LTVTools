function m = subsref(m,L)

switch L(1).type
    case '.'
        L1s = L(1).subs;
        if any(strcmpi(L1s,{'Data','Time','InterpolationMethod'...
                'isTimeInvariant','SplineData','Ts','TimeUnit','StateName','StateUnit',...
                    'InputName','InputUnit','InputGroup','OutputName',...
                    'OutputUnit','OutputGroup','Notes','UserData','Name'}))
            m = m.(L1s);
        elseif any(strcmpi(L1s,{'A','B','C','D'}))
            % Return state matrices as TVMAT
            if m.isTimeInvariant
                m = tvmat( m.Data.(L1s) );
            else
                if isct(m)
                    m = tvmat(m.Data.(L1s),m.Time,m.InterpolationMethod);
                else
                    m = tvmat(m.Data.(L1s),m.Time,m.Ts);
                end
            end
        else
            try
                % Attempt a subsref on SS data
                m = m.Data.(L1s);
            catch
                error(['No property of the class "tvmat" matches the ' ...
                    'identifier "' L(1).subs '". Use PROPERTIES to get ' ...
                    'the list of properties for this class.']);
            end
        end
    case '()'
        % Handle ()-subsref using general UNOP call
        m = tvunop(@subsref,m,L(1));
        
    otherwise
        error(['{}-like reference is not supported for ' class(m) 'objects.'])
end

if length(L)>1
    m = subsref(m,L(2:end));
end
