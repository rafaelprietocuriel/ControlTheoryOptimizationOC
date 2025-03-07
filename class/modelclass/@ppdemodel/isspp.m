function [b diff]=isspp(ppdeObj,varargin)
%
% ISSPP true if limit set is of saddle type

n=statenum(ppdeObj);
[numseigval,numueigval,numceigval]=characteristics(varargin{:});
diff=[numseigval{:}]-[numueigval{:}]-[numceigval{:}];
b=diff==0;
