function [u,v,w] = UnitBox(u,v,w,varargin)
if ~isempty(u), u = max(u,0);u = min(u,1);end
if ~isempty(v), v = max(v,0);v = min(v,1);end
if ~isempty(w), w = max(w,0);w = min(w,1);end