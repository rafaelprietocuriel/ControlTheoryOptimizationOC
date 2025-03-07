function [pn,pnr]=genpattern4name(ocObj,fn)
%
% GENPATTERN4NAME generates the pattern for the saving string.
%
% Example: 'justcoallogsimple_alpha1_0.75_alpha4_1.1_beta_0.5_theta_0.5_thetat_0.5_P_1_PE_0.3_PC_1_phi0_0.01_phiN_1_c3_30_c2_40_wE0_0_wC0_0_HC0_0'
%   m=stdocmodel('justcoallogsimple');
%   [pn,pnr]=genpattern4name(m,'justcoallogsimple_alpha1_0.75_alpha4_1.1_beta_0.5_theta_0.5_thetat_0.5_P_1_PE_0.3_PC_1_phi0_0.01_phiN_1_c3_30_c2_40_wE0_0_wC0_0_HC0_0');
%   Output:
%   pn='alpha1,alpha4,beta,theta,thetat,P,PE,PC,phi0,phiN,c3,c2,wE0,wC0,HC0'
%   pnr='alpha1,alpha4,beta,theta,thetat,P,PE,PC,phi0,phiN,c3,c2,wE0,wC0'
pn='';
pnr='';

o=regexp(fn,'_','split');
if ~strcmp(o{1},modelname(ocObj))
    return
end
parn=parametername(ocObj);
o([1 3:2:length(o)])=[];

rmidx=zeros(1,length(o));
for ii=1:length(o)
    if ~any(strcmp(o{ii},parn))
        rmidx(ii)=ii;
    end
end
pn=cell2vectorstring(o);
pn=pn(2:end-1);
o(rmidx~=0)=[];
pnr=cell2vectorstring(o);
pnr=pnr(2:end-1);
