function b=isfocus(varargin)
%

[numseigval,numueigval,numceigval,seigenvector,seigenvalue,ueigenvector,ueigenvalue,ceigenvector,ceigenvalue]=characteristics(varargin{:});

b=zeros(1,length(seigenvalue));
for ii=1:length(seigenvalue)
    b(ii)=any(abs(imag(seigenvalue{ii}))>0);
end
