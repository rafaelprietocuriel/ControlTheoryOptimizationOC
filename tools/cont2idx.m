function idx=cont2idx(contpar,parval)
%
% CONT2IDX returns indices where a curve crosses a specific parameter value
%
% IDX=CONT2IDX(CONTPAR,PARVAL) returns the index/indices where the
% curve stored as a vector of discrete numbers CONTPAR, usually coming from
% a continuation process, crosses the parameter value PARVAL.

idx=[];
switch numel(parval)
    case 1
        if length(contpar)<2
            return
        end

        if isempty(parval)
            return
        end
        if any(abs(imag(contpar))>0)
            warning('Parameter values have imaginary parts. ')
            contpar=real(contpar);
        end
        [dum,dum,idx]=intersections(1:length(contpar),contpar,1:length(contpar),repmat(parval,1,length(contpar)));
        idx=round(idx);
        if isempty(idx)
            [dum,idx]=min(abs(contpar-parval));
        end
        %         idx0=find(sign(contpar-parval)==0);
        %         idx=find(diff(sign(contpar-parval)));
        %         if ~isempty(idx)
        %             idx=idx+1;
        %         end
        %         idx=unique([idx0 idx]);
    otherwise
        veccontpar=contpar-parval(:,ones(1,size(contpar,2)));
        veccontpar(:,end)=[];
        veccont=diff(contpar,[],2);
        normvec=sqrt(sum(veccont.^2));
        veccont=veccont./normvec(ones(numel(parval),1),:);
        idx=cont2idx(sum(veccont.*veccontpar),0);
end