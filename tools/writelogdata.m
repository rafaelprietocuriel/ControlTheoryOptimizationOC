function logstr=writelogdata(showvariable,existidx)

logstr='';

for ii=1:length(showvariable)
    if ii<length(showvariable)
        eos='\t';
    else
        eos='\n';
    end
    if existidx(ii)
        val=evalin('caller',showvariable{ii});
        if numel(val)>1
            val=max(val(:));
        end
        switch class(val)
            case 'char'
                logstr=[logstr ,[sprintf('%s',val) eos]];
            case {'int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
                logstr=[logstr ,[sprintf('%d',val) eos]];
            case {'double','double'}
                if rem(val,1)==0
                    logstr=[logstr ,[sprintf('%d',val) eos]];
                elseif abs(showvariable{ii})<1e-4
                    logstr=[logstr ,[sprintf('%f',val) eos]];
                else
                    logstr=[logstr ,[sprintf('%e',val) eos]];
                end
        end
        existidx(ii)=1;
    else
        logstr=[logstr ,['[]' eos]];
    end
end
