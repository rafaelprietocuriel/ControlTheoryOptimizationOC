function writeexplicitcontroldynamicsjacobian(m,arcid,filename)

ocStruct=loadmodeldata(modelname(m));
varStruct=generatevariablestruct4standardmodel(ocStruct,[],'INFODETAILS');
varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'PARAMETERVALUES');
varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'IMPLICITINDEX');
fn=fieldnames(varStruct);

filename=[modelname(m) filename];
fid=fopen([filename '.m'],'w');
fprintf(fid,'function [J,coord]=%s(t,depvar,par,arcid)\n',filename);


for ii=1:length(fn)
    if ~varStruct.(fn{ii}).arcdependent
        for jj=1:length(varStruct.(fn{ii}).string)
            fprintf(fid,'%s\n',varStruct.(fn{ii}).string{jj});
        end
    end
end
fprintf(fid,'\nJ=[];\n');
fprintf(fid,'ctrl=%s(t,depvar,par,arcid);\n',[modelname(m) 'OptimalControl']);
fprintf(fid,'lagmcc=%s(t,depvar,par,arcid);\n\n',[modelname(m) 'LagrangeMultiplier']);
fprintf(fid,'switch arcid\n');
for ii=arcid;
    [dUdX,dUdtsym,JU]=generateexplicitcontroldynamics(m,ii);
    transterm=term2standardvariable(m,JU,0,1,0);
    fprintf(fid,'\tcase %d\n',ii);
    fprintf(fid,'\t\tJ=[');
    for jj=1:length(transterm);
        if jj==1
            fprintf(fid,'%s',transterm{jj}(2:end-1));
        else
            fprintf(fid,'\t\t\t%s',transterm{jj}(2:end-1));
        end
        if jj<length(transterm)
            fprintf(fid,'; ...\n');
        else
            fprintf(fid,'];\n');
        end
    end
    fprintf(fid,'\t\tcoord=%s;\n',varStruct.IMPLICITINDEX.string{ii+1});
end
fprintf(fid,'end\n\n');
%fprintf(fid,'J=J(coord,:);\n');
%fprintf(fid,'J=J(:,[%s %s+coord]);',statecostatecoord,statecostatenum);
fclose(fid)