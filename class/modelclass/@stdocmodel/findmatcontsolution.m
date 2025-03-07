function dynPrim=findmatcontsolution(ocObj,biftype,varargin)
%
% 
dynPrim=[];
if isempty(ocObj)
    return
end
matRes=matcontresult(ocObj);
biftype=upper(biftype);
counter=0;
for ii=1:length(matRes)
    for jj=1:length(matRes{ii}.ContinuationInformation)
        b=0;
        if strcmp(deblank(matRes{ii}.ContinuationInformation(jj).label),biftype)
            switch biftype
                case 'LP'
                    if strcmp(matRes{ii}.ContinuationInformation(jj).msg,'Limit point')
                        b=1;
                    else
                        b=0;
                    end
                case 'H'
                    if strcmp(matRes{ii}.ContinuationInformation(jj).msg,'Hopf')
                        b=1;
                    else
                        b=0;
                    end
            end
            if b
                counter=counter+1;
                idx=matRes{ii}.ContinuationInformation(jj).index;
                par=matRes{ii}.ContinuationSolution.modelparameter;
                aparidx=matRes{ii}.ContinuationSolution.userinfo.varyparameterindex;
                epstruct.arcarg=matRes{ii}.ContinuationSolution.arcarg(idx);
                epstruct.y=matRes{ii}.ContinuationSolution.y(:,idx);
                apar=matRes{ii}.ContinuationSolution.userinfo.varyparametervalue(idx);
                par(aparidx)=apar;
                ocObjN=changeparametervalue(ocObj,par);
                epstruct.modelparameter=par;
                epstruct.userinfo.activeparameter=apar;
                epstruct.userinfo.activeparameterindex=aparidx;
                epstruct.userinfo.stdocmodel=ocObjN;
                dynPrim{counter}=dynprimitive(epstruct,ocObjN);
                dynPrim{counter}.linearization=linearize(dynPrim{counter},ocObjN,'dependentvar',1);
            end
        end
    end
end
