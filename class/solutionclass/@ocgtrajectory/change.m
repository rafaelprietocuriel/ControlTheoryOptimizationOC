function ocgTrj=change(ocgTrj,varargin)

for ii=2:2:nargin
    try
        ocgTrj.(varargin{ii-1})=varargin{ii};
    catch
        if strcmp(varargin{ii-1},'arcarg')
            oldarcarg=ocgTrj.octrajectory.(varargin{ii-1});
            ctrl=control(ocgTrj);
            x=state(ocgTrj);
            l=costate(ocgTrj);
            y=ocgTrj.octrajectory.y;
            ocObj=loadmodel(ocgTrj);
            diffidx=find(oldarcarg-varargin{ii});
            for jj=1:length(diffidx)
                actargold=oldarcarg(diffidx(jj));
                actargnew=varargin{ii}(diffidx(jj));
                olddim=canonicalsystemdimension(ocObj,actargold);
                newdim=canonicalsystemdimension(ocObj,actargnew);
                coord=implicitcontrolcoordinate(ocObj,actargnew);
                yycs=y{diffidx(jj)}(1:olddim,:);
                yycsa=y{diffidx(jj)}(olddim+1:end,:);
                yycsn=zeros(newdim,size(yycs,2));
                yycsn(1:newdim,:)=[x{diffidx(jj)};l{diffidx(jj)};ctrl{diffidx(jj)}(coord,:)];
                ocgTrj.octrajectory.y{diffidx(jj)}=[yycsn;yycsa];
                ocgTrj.odenum(diffidx(jj))=newdim;
            end
            ocgTrj.octrajectory.(varargin{ii-1})=varargin{ii};
        elseif strcmp(varargin{ii-1},'state')
            n=statenum(ocgTrj);
            for jj=1:arcnum(ocgTrj)
                ocgTrj.octrajectory.y{jj}(1:n,:)=varargin{ii}{jj};
            end
        elseif strcmp(varargin{ii-1},'modelparameter')
            ocgTrj.octrajectory.modelparameter=varargin{ii};
        end

    end
end
