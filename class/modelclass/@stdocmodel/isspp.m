function [b,diff]=isspp(ocObj,varargin)
%
% ISSPP true if limit set is of saddle type

numseigval=characteristics(varargin{:});
diff=[numseigval{:}]-statenum(ocObj);
b=diff==0;
