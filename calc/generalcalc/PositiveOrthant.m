function [u,v,w] = PositiveOrthant(u,v,w,varargin)
if ~isempty(u), u = max(u,0);end
if ~isempty(v), v = max(v,0);end
if ~isempty(w), w = max(w,0);end