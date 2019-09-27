function a = getSplineData(Data,T)
% getSplineData Compute cubic spline coefficients
%
% A=getSplineData(Data,Time) computes cubic spline coefficients to 
% interpolate the input Data and Time.  Time is an Nt-by-1 vector of
% increasing time values and Data is an 1-by-1-by-Nt double array. The
% spline is given by cubic functions on the intervals [Time(k), Time(k+1)] 
% for k=1,...,Nt-1. Thus the spline is defined by a total of 4*(Nt-1)
% coefficients (Nt-1 intervals and 4 coefficients per interval).
% The cubic functions are defined to interpolate { Time(k), 
% Data(1,1,k) }_{k=1}^Nt and satisfy constraints to ensure continuity of 
% the function and its first/second derivatives.  The output A is an 
% 4*(Nt-1)-by-1 vector of spline coefficients. The function EVALT is used
% to evaluate the spline at any time point.
%
% If Data is N1-by-N2-by-...-by-Nad-by-Nt then A defines spline 
% coefficients as a 4*(Nt-1)-by-(N1*N2*...*Nad) matrix.
%
% A = getSplineData(A) fills the SplineData property of a TVMAT A.
% This allows for more efficient interpolation by EVALT.

% The notation in this file is based on the following reference:
%   Moore, "Finite Horizon Robustness Analysis using IQCs," MS, 2015.


%% Handle TVMAT Input
% SplineData is initialized to be empty and it is computed when calling
% EVALT (if InterpolationMethod = 'Spline').  The code below allows
% SplineData to be pre-computed once for efficient repeated calls to EVALT.
if nargin==1 && isa(Data,'tvmat')
    a = Data;
    a.SplineData = getSplineData(a.Data,a.Time);
    return
end


%% Form Linear Equations
% The cubic spline in the k^th interval [T(k),T(k+1)] is 
%    fk(t) = (a0)_k + (a1)_k*dt + (a2)_k*dt^2 + (a3)_k dt^3
% where dt = t-Tk.  The constraints in the reference yield:
%    (a0)_k = pk   for k=1,...,N-1
%     H*x = b 
% where x = [(a1)_1; (a2)_1; (a3)_1; ....; 
%            (a1)_{N-1}; (a2)_{N-1}; (a3)_{N-1}];
% and pk is the k^th value in Data.
Nt = numel(T);    % Number of interpolation points
H = zeros(3*(Nt-2)+1,3*(Nt-1));
for k=1:Nt-2
    hk = T(k+1)-T(k);
    Hk = [hk hk^2 hk^3 0 0 0; 1 2*hk 3*hk^2 -1 0 0; 0 2 6*hk 0 -2 0];    
    H( 3*(k-1)+(1:3), 3*(k-1)+(1:6)) =  Hk;
end

hk = T(Nt)-T(Nt-1);
Hk = [0 0 0 hk hk^2 hk^3];
H(end,end-5:end) = Hk;
%H = [H; zeros(1,3*(Nt-3)) Hk];

if Nt==2
    % Natural Conditions
    H = [hk hk^2 hk^3; 0 2 0;0 2 6*hk];
elseif Nt==3
    % Natural Conditions
    Hk = [0 2 0 0 0 0;0 0 0 0 2 6*hk];
    H = [H; Hk zeros(1,3*(Nt-3))];    
else
    % Not-a-knot (for N>=4)    
    Hk = [0 0 6 0 0 -6];
    H = [H; zeros(1,3*(Nt-3)) Hk];
    Hk = [0 0 6 0 0 -6];
    H = [H; Hk zeros(1,3*(Nt-3))];    
end

%% Generate Spline Basis Coefficients

% Data can be an N1-by-N2-...-by-Na-by-Nt array. Reshape and transpose
% to represent as a two-dimensional Nt-by-Ne array.
sz = size(Data);
nad = numel(sz)-3;
id = repmat({':'},1,nad);
Ne = numel( Data(:,:,id{:},1) );
p = reshape(Data,[Ne,Nt])';

% Solution vector for least squares problem
b = zeros(3*(Nt-1),Ne);
for k=1:Nt-1
    idx = 3*(k-1)+(1:3);
    b(idx,:) = [p(k+1,:)-p(k,:); zeros(2,Ne)];
end

% Solve least squares
% x is 3*(Nt-1)-by-Ne and contains all coefficients except the constants.
x = H\b;

% Incorporate constants into least squares solution
% a is 4*(Nt-1)-by-Ne and contains all coeffs including the constants.
a = zeros(4*(Nt-1),Ne);
for k=1:Nt-1
    xidx = 3*(k-1)+(1:3);
    aidx = 4*(k-1)+(1:4);
    a( aidx, : ) = [p(k,:); x(xidx,:)];
end

