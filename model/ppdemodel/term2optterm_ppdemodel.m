function optterm=term2optterm_ppdemodel(ocStruct,term,arcid,opt)
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
ctrlnum=retrieveppdemodelinformation(ocStruct,'controlnum');

inequalitycontrolconstraint=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraint');
ccflag=~isempty(inequalitycontrolconstraint.value);
controlname=retrieveppdemodelinformation(ocStruct,'controlname');
lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
if ~iscell(term)
    term=cell(term);
end
optterm=cell(length(term),1);

%maximizingsolution=retrieveppdemodelinformation(ocStruct,'maximizingsolution',arcid,symkernel);
ctrl=ocStruct.foc.value.control.(arcidentifier2field(arcid)).term;
if ccflag
    lmcc=ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcid)).term;
end
for ii=1:length(term)
    optterm{ii}=term{ii};
    for jj=1:ctrlnum.value
        optterm{ii}=ocmatsubs(optterm{ii},[controlname.value{jj} '=' ctrl{jj}],symkernel);
        if Simplify
            optterm{ii}=ocmatsimple(optterm{ii},symkernel);
        end
    end
    if ccflag
        for jj=1:length(lagrangemultipliercontrolname.value)
            optterm{ii}=ocmatsubs(optterm{ii},[lagrangemultipliercontrolname.value{jj} '=' lmcc{jj}],symkernel);
            if Simplify
                optterm{ii}=ocmatsimple(optterm{ii},symkernel);
            end
        end
    end
end

