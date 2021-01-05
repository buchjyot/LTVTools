function m = subsref(m,L)

switch L(1).type
    case '.'
        L1s = L(1).subs;
        if any(strcmpi(L1s,{'Data','Ts','Time','TimeUnit','InterpolationMethod',...
                'isTimeInvariant','SplineData'}))
            m = m.(L1s);
        elseif any(strcmpi(L1s,{'A','B','C','D'}))
            % Return state matrices as TVMAT
            if m.isTimeInvariant
                m = tvmat( m.Data.(L1s) );
            else
                m = tvmat(m.Data.(L1s),m.Time,m.InterpolationMethod);
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
