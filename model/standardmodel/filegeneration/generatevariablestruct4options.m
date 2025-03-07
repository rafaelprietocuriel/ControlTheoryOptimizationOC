function varStruct=generatevariablestruct4options(opt,varargin)
%
jacobianflag=getocoptions(opt,'INIT','Jacobian');
varnamelist={'JACOBIANNUMERICALFLAG','JACOBIANEXPLICITFLAG','JACOBIANIMPLICITFLAG'};

for ii=1:length(varnamelist)
    varname=varnamelist{ii};
    varStruct.(varname).arcdependent=0;
    varStruct.(varname).arcidentifier=[];
    varStruct.(varname).vectorize=0;
    varStruct.(varname).type='';
    varStruct.(varname).multline=0;

    switch varname

        case 'JACOBIANNUMERICALFLAG'
            varStruct.(varname).string=num2str(strcmp(jacobianflag,'numerical'));

        case 'JACOBIANEXPLICITFLAG'
            varStruct.(varname).string=num2str(strcmp(jacobianflag,'explicit'));

        case 'JACOBIANIMPLICITFLAG'
            varStruct.(varname).string=num2str(strcmp(jacobianflag,'implicit'));
    end
end