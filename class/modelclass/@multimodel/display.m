function display(mmObj)

if isequal(get(0,'FormatSpacing'),'compact')
    disp([inputname(1) ' =']);
    disp('ocmatclass: multimodel')
    disp(['    modelname : ' modelname(mmObj)]);
    disp('consisting of ')
else
    disp([inputname(1) ' =']);
    disp('ocmatclass: multimodel')
    disp(' ');
    disp(['    modelname : ' modelname(mmObj)]);
    disp(' ')
    disp('consisting of ')
end
for ii=1:numberofmodels(mmObj)
    display(mmObj.Model{ii})
%     if isequal(get(0,'FormatSpacing'),'compact')
%         disp([inputname(1) '(' num2str(ii) ') =']);
%         if ~isempty(mmObj(ii))
%             [parval,parvar]=parametervalue(mmObj.Model{ii});
%             disp('ocmatclass: stdocmodel')
%             disp(['    modelname : ' modelname(mmObj.Model{ii})]);
%             for jj=1:parameternum(mmObj.Model{ii})
%                 disp(['    ' parvar{jj} ' : ' num2str(parval(jj))]);
%             end
%         else
%             disp('empty stdocmodel');
%         end
%     else
%         disp(' ')
%         disp([inputname(1) '(' num2str(ii) ') =']);
%         disp(' ')
%         if ~isempty(mmObj.Model{ii})
%             [parval,parvar]=parametervalue(mmObj.Model{ii});
%             disp('ocmatclass: stdocmodel')
%             disp(' ');
%             disp(['    modelname : ' modelname(mmObj.Model{ii})]);
%             disp(' ')
%             for jj=1:parameternum(mmObj.Model{ii})
%                 disp(['    ' parvar{jj} ' : ' num2str(parval(jj))]);
%             end
%         else
%             disp('empty stdocmodel');
%         end
%         disp(' ');
%     end
end