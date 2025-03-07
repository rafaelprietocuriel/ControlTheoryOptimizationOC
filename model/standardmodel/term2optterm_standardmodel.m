function optterm=term2optterm_standardmodel(ocStruct,term,arcid,opt)
%
%
if isempty(ocStruct)
    optterm=term;
    return
end

if nargin==3
    opt=defaultocoptions;
end


Simplify=strcmp(getocoptions(opt,'INIT','Simplify'),'on');
symkernel=getsymkernel;

if ~ischar(arcid)
    arcid=num2str(arcid);
end
tmp=retrievemodelinformation(ocStruct,'controlnum');
ctrlnum=tmp.value;

tmp=retrievemodelinformation(ocStruct,'inequalitycontrolconstraint');
inequalitycontrolconstraint=tmp.value;
ccflag=~isempty(inequalitycontrolconstraint);
tmp=retrievemodelinformation(ocStruct,'controlname');
controlname=tmp.value;
tmp=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
lagrangemultipliercontrolname=tmp.value;
if ~iscell(term)
    term=cell(term);
end
optterm=cell(length(term),1);

%maximizingsolution=retrievemodelinformation(ocStruct,'maximizingsolution',arcid,symkernel);
ctrl=ocStruct.foc.value.control.(arcidentifier2field(arcid)).term;
if ccflag
    lmcc=ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcid)).term;
end
for ii=1:length(term)
    optterm{ii}=term{ii};
    for jj=1:ctrlnum
        optterm{ii}=ocmatsubs(optterm{ii},[controlname{jj} '=' ctrl{jj}],symkernel);
        if Simplify
            optterm{ii}=ocmatsimple(optterm{ii},symkernel);
        end
    end
    if ccflag
        for jj=1:length(lagrangemultipliercontrolname)
            optterm{ii}=ocmatsubs(optterm{ii},[lagrangemultipliercontrolname{jj} '=' lmcc{jj}],symkernel);
            if Simplify
                optterm{ii}=ocmatsimple(optterm{ii},symkernel);
            end
        end
    end
end

