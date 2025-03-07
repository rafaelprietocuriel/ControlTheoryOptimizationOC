function out=octype(ocgMTrj)
% octype returns the type of the ocgradmultipath, possible types are: 
% concentrated ... non-distributed variables
% age ... age distributed
% space ... spatial distributed

out=ocgMTrj.type;