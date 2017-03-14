%function that finds the coordinates of a node in a cartesian grid.
%Coordinates do not depend in the size of the grid in z, because we assume
%that the point is not out of the grid.
%Inputs are the id of the node (in Matlab notation: y, then x, then z),
%then the size of the grid in y and x.

function [y,x,z] = findcoord(id,sizey,sizex)

z = ceil(id./(sizey.*sizex));
id2d = (id-(sizey.*sizex.*(z-1)));
x = ceil(id2d./sizey);
y = id2d-((x-1).*sizey);