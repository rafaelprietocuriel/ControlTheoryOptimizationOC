function writeexplicitcontroldynamics(m,arcid,filename)

ocStruct=loadmodeldata(modelname(m));
varStruct=generatevariablestruct4standardmodel(ocStruct,[],'INFODETAILS');
varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'PARAMETERVALUES');
varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'IMPLICITINDEX');
fn=fieldnames(varStruct);

filename=[modelname(m) filename];
fid=fopen([filename '.m'],'w');
fprintf(fid,'function [dudt,coord]=%s(t,depvar,par,arcid)\n',filename);


for ii=1:length(fn)
    if ~varStruct.(fn{ii}).arcdependent
        for jj=1:length(varStruct.(fn{ii}).string)
            fprintf(fid,'%s\n',varStruct.(fn{ii}).string{jj});
        end
    end
end
fprintf(fid,'\ndudt=[];\ncoord=[];\n');
fprintf(fid,'ctrl=%s(t,depvar,par,arcid);\n',[modelname(m) 'OptimalControl']);
fprintf(fid,'lagmcc=%s(t,depvar,par,arcid);\n\n',[modelname(m) 'LagrangeMultiplier']);
fprintf(fid,'switch arcid\n');
for ii=arcid;
    [dUdX,dUdtsym]=generateexplicitcontroldynamics(m,ii);
    transterm=term2standardvariable(m,dUdtsym,0,1,1);
    fprintf(fid,'\tcase %d\n',ii);
    fprintf(fid,'\t\tdudt=[');
    for jj=1:length(transterm);
        if jj==1
            fprintf(fid,'%s',transterm{jj});
        else
            fprintf(fid,'\t\t\t%s',transterm{jj});
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
fprintf(fid,'dudt=dudt(coord,:);\n');
fclose(fid)