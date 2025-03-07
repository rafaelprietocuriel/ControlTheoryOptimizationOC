function out=spacegeometrytype(ppdeObj,varargin)

out=[];
if isempty(ppdeObj)
    return
end

if isfield(ppdeObj.Model.spacegeometry,'R1') && isfield(ppdeObj.Model.spacegeometry.R1,'interval')
    out='interval';
end