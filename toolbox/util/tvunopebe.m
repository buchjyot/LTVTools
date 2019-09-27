function B = tvunopebe(opfh,A,varargin)
%% TVUNOPEBE Unary operations applied element-by-element for a TVMAT
%
% TVUNOP is a utility function that handles unary operations (e.g. transpose,
% ctranspose, etc.) on TVMATs. The user should not call TVUNOP directly.
% Instead each unary operation has a front end function that the user
% calls (e.g. transpose). The front end function then call calls
% TVUNOP  to do the required calculations/operations.

% Perform unary operation at each point 
B = A;
B.Data = opfh(A.Data,varargin{:});

