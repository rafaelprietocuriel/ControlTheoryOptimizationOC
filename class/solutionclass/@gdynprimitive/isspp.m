function [b,diff]=isspp(gdynPrim)
%
% ISSPP true if limit set is of saddle type

ocObj=loadmodel(gdynPrim);
numseigval=characteristics(gdynPrim);
diff=[numseigval{1}]-statenum(ocObj);
b=diff==0;
