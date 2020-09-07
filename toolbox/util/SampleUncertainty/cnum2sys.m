function [Delta,beta] = cnum2sys(delta0,w0)
% Delta = cnum2sys(delta0,w0)
%
% The complex number delta0 and frequency w0>0 are given.
% This function constructs a stable, LTI system Delta(s) such that
% Delta(j w0) = delta0 and || Delta ||_infty <= |delta0|

if isreal(delta0)
    beta = 0;
    Delta = ss([],[],[],delta0);
elseif w0>0 && ~isinf(w0)
    % Polar form with sigma*c*exp(j*phi) where sigma=+/- 1.
    % angle returns phase in rads in the range [-pi,pi]. If 
    % Imag(delta0)>0 then phase is in (0,pi) and if Imag(delta0)<0
    % then phase is in (-pi,0).
    c = abs(delta0);
    phi = angle(delta0);
    sigma = +1;
    if imag(delta0)<0
        % Map phase to [0, phi] and change sign coefficient
        phi = phi + pi;
        sigma = -1;
    end
    
    % Compute pole location for Delta(s)
    beta = w0*tan(phi/2);
    
    % Express Delta in state-space form:
    % Delta(s) = sigma*c*(s-beta)/(s+beta) 
    %          = sigma*c*[1 - (2*beta)/(s+beta)]
    %          = sigma*c - sigma*(2*beta*c)/(s+beta)
    A = -beta;
    B = sqrt(2*beta*c);
    C = -sigma*B;
    D = sigma*c;
    Delta = ss(A,B,C,D);
else
    error('If delta0 is not real then w0 must be >0 and finite');
end
