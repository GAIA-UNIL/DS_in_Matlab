function coord = find_node_in_grid(points , grid , no_nodes)
for i=1:no_nodes
    coord(i,1)=grid(points(i,1),points(i,2));
end