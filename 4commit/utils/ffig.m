function [h]=ffig(position)
%plot a giant figure. 70% of screen height/width by default. 
%  Assumes position is a vector of the following format.
%	[ xpos ypos xsize ysize ], 
%  where xpos/ypos is the x and y distance of the top left corner of 
%  the figure from the top left corner of the screen, and xsize/ysize 
%  is the proportion of the screen size (0->1) that you want the figure
%  to cover.

if ~exist('position','var')||isempty(position);
    position=[ .15 .15 .7 .7 ];
end

h=figure('units','normalized','position',position);
end
