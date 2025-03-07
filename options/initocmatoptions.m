function out = initocmatoptions

% OCMAT Toolbox options
%
% Field AutoDataFileGeneration: y/n/(a) (yes / no / ask (default))
% Field OutputFormatODE 'ocmat' (default) / 'matlab'

mainoptionsdir=getocmatpath('options');

try
    load OCMatOptions
    out=ocmatoptions;
catch

    ocmatoptions.AutoDataFileGeneration='a';
    ocmatoptions.OutputFormatODE='ocmat';
    out=ocmatoptions;
    save(fullfile(mainoptionsdir,['OCMatOptions.mat']),'ocmatoptions')
    
end