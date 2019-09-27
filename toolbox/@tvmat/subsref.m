function m = subsref(m,L)

switch L(1).type
    case '.'
        L1s = L(1).subs;        
        if any(strcmpi(L1s,{'Data','Time','InterpolationMethod',...
                'isTimeInvariant','SplineData','Ts','TimeUnit','Name'}))
            m = m.(L1s);
        else
            error(['No property of the class "tvmat" matches the ' ...
                'identifier "' L(1).subs '". Use PROPERTIES to get ' ...
                'the list of properties for this class.']); 
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
