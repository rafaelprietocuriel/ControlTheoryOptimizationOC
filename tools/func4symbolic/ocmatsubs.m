function NEWf=ocmatsubs(OLDf,expr,symkernel)
%ocmatsubs
NEWf=OLDf;
if nargin==2
    symkernel='maple';
end

switch symkernel
    case 'maple'
        NEWf=maple('subs',expr,OLDf);
    case 'mupad'
        expr=regexp(removeouterbrackets(expr),',','split');
        for ii=1:numel(expr)
            splitexpr=regexp(expr{ii},'=','split');
            LS{ii}=splitexpr{1};
            RS{ii}=splitexpr{2};
        end
        %NEWf=char(simple(subs(sym(OLDf),{LS{:}},{RS{:}},0)));
        if verLessThan('symbolic','8')
            try
                NEWf=subs(sym(OLDf),{LS{:}},{RS{:}},0);
            catch
                NEWf=subs(sym(OLDf),{LS{:}},{RS{:}});
            end
            if isnumeric(NEWf)
                NEWf=num2str(NEWf);
            else
                NEWf=char(NEWf);
            end
        else
            LS=cellfun(@mystr2sym,LS,'UniformOutput',false);
            RS=cellfun(@mystr2sym,RS,'UniformOutput',false);
            ischar=false;
            if ~isa(OLDf,'sym')
                OLDf=mystr2sym(OLDf);
                ischar=true;
            end
            try
                NEWf=subs(OLDf,{LS{:}},{RS{:}},0);
            catch
                NEWf=subs(OLDf,{LS{:}},{RS{:}});
            end
            if isnumeric(NEWf)
                NEWf=num2str(NEWf);
            else
                if ischar
                    NEWf=char(NEWf);
                end
            end
        end
    otherwise
        ocmatmsg('No kernel for the symbolic toolbox specified.\n')
end

