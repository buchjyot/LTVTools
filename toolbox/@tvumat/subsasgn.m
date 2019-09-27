function m = subsasgn(m,L,RHS)

switch L(1).type
    case '.'
        % A.B.C = RHS  -->  TMP = A.B; TMP.C = RHS; A.B = TMP:
        % XXX Code below does not work on "Function" datatypes
        try
            if length(L) == 1
                tmp = RHS;
            else
                tmp = subsref(m,L(1));
                tmp = subsasgn(tmp,L(2:end),RHS);
            end
            L1s = L(1).subs;
            if strcmpi(L1s,'Time')
                % Store Time as a column vector
                m.Time = tmp(:);
                
                % Reset spline data
                m.SplineData = [];
            elseif strcmpi(L1s,'Data')
                m.Data = tmp;
                
                % Reset spline data
                m.SplineData = [];
            elseif strcmpi(L1s,'InterpolationMethod')
                m.InterpolationMethod = tmp;
                
                % Reset spline data
                m.SplineData = [];
            elseif strcmpi(L1s,'SplineData')
                % XXX This allows access to SplineData so that it can
                % be set by getSplineData in LTVUtil. It is probably
                % not good to make this publicly accessible.
                m.SplineData = tmp;
            elseif any(strcmpi(L1s,{'A','B','C','D'}))
                % Assign underlying state-space data
                % XXX What if Time/InterpolationMethod is inconsistent?
                if isa(tmp,'TVMAT')
                    tmp = tmp.Data;
                end
                m.Data.(L1s) = tmp;
            else
                try
                    % Attempt a subsasgn on SS data
                    m.Data.(L1s) = tmp;
                catch
                    error(['No property of the class "tvmat" matches the ' ...
                        'identifier "' L1s '". Use PROPERTIES to get ' ...
                        'the list of properties for this class.']);
                end
            end
        catch ME
            rethrow(ME);
        end
    case '()'
        % Handle ()-subsasgn using general BINOP call
        L1 = L(1);
        opfh = @(A,B) subsasgn(A,L1,B);
        m = tvbinop(opfh,m,RHS);
    case '{}'
        error(['{}-like {} SUBSASGN not supported for ' class(m) 'objects.'])
end

% Check validity
[pflag,errstr] = isvalid(m);
if pflag==0
    error(errstr);
end

