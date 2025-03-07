function [idx,solObj]=findimpulsesolution(ocObj,numintimpulse,numinitimpulse,numendimpulse,numfiximpulse,varargin)
%
% 

solObj=[];
if nargin==2
    numinitimpulse=0;
    numendimpulse=0;
    numfiximpulse=0;
end
if nargin==3
    numendimpulse=0;
    numfiximpulse=0;
end
if nargin==4
    numfiximpulse=0;
end
if isempty(ocObj)
    return
end
ocEx=extremalsolution(ocObj);

idx=zeros(1,length(ocEx));
counter=0;
for ii=1:length(ocEx)
    jumparg=jumpargument(ocEx{ii});
    switch class(ocEx{ii})
        case 'hybridoctrajectory'
            s1=sum(jumparg==1);
            s2=(jumparg(1)==-1);
            s3=(jumparg(end)==-1);
            s4=sum(jumparg(2:end-1)==-1);
            if s1==numintimpulse && s2==numinitimpulse && s3==numendimpulse && s4==numfiximpulse
                counter=counter+1;
                idx(counter)=ii;
            end
    end
end
idx(counter+1:end)=[];
solObj=ocEx(idx);