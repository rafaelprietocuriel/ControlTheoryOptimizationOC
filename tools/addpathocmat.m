function addpathocmat

mainpath=getocmatpath('numtools');
mainpath=fullfile(mainpath,'ocmat');

addpath(genpath(mainpath))
rmpath(genpath(fullfile(mainpath,'doc')))
rmpath(genpath(fullfile(mainpath,'general','adaptedmatcontfiles')))
rmpath(genpath(fullfile(mainpath,'general','dasslc')))
rmpath(genpath(fullfile(mainpath,'general','impulse')))
rmpath(genpath(fullfile(mainpath,'general','manifolds')))
rmpath(genpath(fullfile(mainpath,'model','out')))
rmpath(genpath(fullfile(mainpath,'model','template')))
rmpath(genpath(fullfile(mainpath,'_gsdata_')))
savepath 