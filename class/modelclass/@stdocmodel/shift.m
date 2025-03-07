function ocEPs=shift(ocObj,k,varargin)


if nargin==1
    ocEPs=dynprimitive([]);
    return
end
numx=statenum(ocObj);
N=parametervalue(ocObj,'N');
if numx>N+1
    ocEPs=dynprimitive([]);
    return
end
if isempty(k)
    k=ceil(N/2);
end
ocEPs=cell(1,nargin-2);
if ~mod(N,2) && k==ceil(N/2)
    for ii=1:nargin-2
        x=state(ocObj,varargin{ii});
        l=costate(ocObj,varargin{ii});
        l=[l(k+1:end);l(2:k+1)];
        l([k+1])=l([k+1])*2;
        l([1 end])=l([1 end])/2;
        ocEPStruct.octrajectory=struct(varargin{ii}.octrajectory);
        ocEPStruct.octrajectory.y=[x(k+1:end);x(2:k+1);l];
        ocEPStruct.linearization=[];
        ocEPs{ii}=dynprimitive(ocEPStruct);
    end
end