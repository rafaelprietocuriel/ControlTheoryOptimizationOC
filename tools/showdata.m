function varargout=showdata(ocEP)

m=stdocmodel(modelname(ocEP),[],[],0);

par=modelparameter(ocEP);
m=changeparametervalue(m,par);
parn=parametername(m);

staten=state(m);
if isa(ocEP,'gdynprimitive')
    lm=lagrangemultiplier(ocEP);
    lm=lm{1};
    cstr=constraint(ocEP);
    cstr=cstr{1};

    ctrln=string2cell(removematrixstring(char(control(m))),'vector');
    x=state(ocEP);
    x=x{1};

    ctrl=control(ocEP);
    ctrl=ctrl{1};
else
    lm=lagrangemultiplier(m,ocEP);
    cstr=constraint(m,ocEP);
    ctrln=string2cell(removematrixstring(char(control(m))),'vector');
    x=state(m,ocEP);

    ctrl=control(m,ocEP);
end
arcid=arcargument(ocEP);
ccomb=constraintcombination(m,arcid);
id=constraintidentifier(m);
cstreq=constraint(m);
lmcc=lagrangemultiplier(m);
fprintf('\nFormal Characteristics:\n');
fprintf('\nArcidentifier : %d\n',arcid);
fprintf('Active constraint(s) : ');
for ii=1:length(ccomb)
    if ii<length(ccomb)
        fprintf('%s,',ccomb{ii});
    else
        fprintf('%s\n',ccomb{ii});
    end
end
fprintf('\nStates:\n');
for ii=1:length(staten)
    if mod(x(ii),1)==0
        fprintf('%s : %d\n',staten{ii},x(ii));
    else
        fprintf('%s : %f\n',staten{ii},x(ii));
    end
end

fprintf('\nControls:\n');
for ii=1:length(ctrln)
    if mod(ctrl(ii),1)==0
        fprintf('%s : %d\n',ctrln{ii},ctrl(ii));
    else
        fprintf('%s : %f\n',ctrln{ii},ctrl(ii));
    end
end

fprintf('\nConstraints:\n');
for ii=1:length(cstr)
    if mod(cstr(ii),1)==0
        fprintf('%s (%s>=0) : %d\n',id{ii},char(cstreq(ii)),cstr(ii));
    else
        fprintf('%s (%s>=0) : %f\n',id{ii},char(cstreq(ii)),cstr(ii));
    end
end
fprintf('\nLagrange Multiplier:\n');
for ii=1:length(cstr)
    if mod(lm(ii),1)==0
        fprintf('%s : %d\n',lmcc{ii},lm(ii));
    else
        fprintf('%s : %f\n',lmcc{ii},lm(ii));
    end
end


fprintf('\nParameter Values:\n');
for ii=1:length(parn)
    if mod(par(ii),1)==0
        fprintf('%s : %d\n',parn{ii},par(ii));
    else
        fprintf('%s : %f\n',parn{ii},par(ii));
    end

end

if nargout>=1
    varargout{1}=cstr;
end
if nargout>=2
    varargout{2}=lm;
end