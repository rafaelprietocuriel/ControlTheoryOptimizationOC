function [u,v,w] = ValueBoundBox(u,v,w,options,varargin)
if ~isempty(u) 
    Siz = size(u);
    if ~isfield(options,'LowerValueU'),options.LowerValueU = options.LowerBoundU;end
    siz = size(options.LowerBoundU);
    if any(siz > 1) && prod(siz) ~= prod(Siz)% bound is neither scalar nor full matrix
        switch siz(1)% index of the first nonsingleton dimention
            case 1
                for k = 1:siz(2)
                    u(:,k,:) = max(u(:,k,:),options.LowerBoundU(1,k)); 
                    for j = 1:Siz(1)
                        u(j,k,u(j,k,:) == options.LowerBoundU(1,k)) = options.LowerValueU(1,k); 
                    end
                end
            otherwise, error('Nontested option in function ValueBoundBox')
        end
    else
        u = max(u,options.LowerBoundU);
        u(u == options.LowerBoundU) = options.LowerValueU; 
    end
    if ~isfield(options,'UpperValueU'),options.UpperValueU = options.UpperBoundU;end
    siz = size(options.UpperBoundU); % UpperBound
    if any(siz > 1) && prod(siz) ~= prod(Siz)% bound is neither scalar nor full matrix
        switch siz(1)% index of the first nonsingleton dimention
            case 1
                for k = 1:siz(2)
                    u(:,k,:) = min(u(:,k,:),options.UpperBoundU(1,k)); 
                    for j = 1:Siz(1)
                        u(j,k,u(j,k,:) == options.UpperBoundU(1,k)) = options.UpperValueU(1,k); 
                    end
                end
            otherwise, error('Nontested option in function ValueBoundBox')
        end
    else
        u = min(u,options.UpperBoundU);
        u(u == options.UpperBoundU) = options.UpperValueU; 
    end
end
if ~isempty(v)
    Siz = size(v);
    if ~isfield(options,'LowerValueV'),options.LowerValueV = options.LowerBoundV;end
    siz = size(options.LowerBoundV); % LowerBound
    if any(siz > 1) && prod(siz) ~= prod(Siz)% index of the first nonsingleton dimention           
                for j = 1:siz(1)
                    v(j,:) = max(v(j,:), options.LowerBoundV(j)); 
                    v(j,v(j,:) == options.LowerBoundV(j)) = options.LowerValueV(j); 
                end
    else
                v = max(v, options.LowerBoundV); 
                if any(siz > 1)
                    v(v == options.LowerBoundV) = options.LowerValueV(v == options.LowerBoundV); 
                else
                    v(v == options.LowerBoundV) = options.LowerBoundV; 
                end
    end
    if ~isfield(options,'UpperValueV'),options.UpperValueV = options.UpperBoundV;end
    siz = size(options.UpperBoundV); % UpperBound
    if any(siz > 1) && prod(siz) ~= prod(Siz)
                for j = 1:siz(1)
                    v(j,:) = min(v(j,:), options.UpperBoundV(j)); 
                    v(j,v(j,:) == options.UpperBoundV(j)) = options.UpperValueV(j); 
                end
    else
                v = min(v, options.UpperBoundV);
                if any(siz > 1)
                    v(v == options.UpperBoundV) = options.UpperValueV(v == options.UpperBoundV); 
                else
                    v(v == options.UpperBoundV) = options.UpperValueV; 
                end
    end
end
if ~isempty(w) 
    Siz = size(w);
    if ~isfield(options,'LowerValueW'),options.LowerValueW = options.LowerBoundW;end
    siz = size(options.LowerBoundW);
    if any(siz ~= 1) && prod(siz) ~= prod(Siz)% bound is neither scalar nor full matrix
        for j = 1:siz(2)
            w(:,j) = max(w(:,j),options.LowerBoundW(1,j)); 
            w(w(:,j) == options.LowerBoundW(j),j) = options.LowerValueW(1,j); 
        end
    else
        w = max(w,options.LowerBoundW);
        w(w == options.LowerBoundW) = options.LowerValueW(w == options.LowerBoundW); 
    end
    if ~isfield(options,'UpperValueW'),options.UpperValueW = options.UpperBoundW;end
    siz = size(options.UpperBoundW); % UpperBound
    if any(siz ~= 1) && prod(siz) ~= prod(Siz)% bound is neither scalar nor full matrix
        for j = 1:siz(2)
            w(:,j) = min(w(:,j,:),options.UpperBoundW(1,j)); 
            w(w(:,j) == options.UpperBoundW(j),j) = options.UpperValueW(1,j); 
        end
    else
        w = min(w,options.UpperBoundW);
        w(w == options.UpperBoundW) = options.UpperValueW(w == options.UpperBoundW); 
    end
end