function varargout = double2tvmat(varargin)
%% DOUBLE2TVMAT
%
% This helper function wraps the incoming data and time into saperate tvmat
%
% Usage: 
% >> [D1,D2] = double2tvmat(Data1,Data2,...,Tgrid)

nin = nargin;
varargout = cell(nin,1);
t = varargin{end};
for i = 1:nin-1
    varargout{i} = tvmat(varargin{i},t);
end