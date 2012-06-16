% a simple example for running pops_to_routers.m
popgraph = struct();% struct containing information about the PoPs.
popgraph.adjacency = [0 1 1 1; 1 0 1 0; 1 1 0 0; 1 0 0 0];% pop adjacency matrix
popgraph.num_access_routers = [2 3 2 2]; % the number of access routers for each pop
rules = struct();% the rules for creating the router level graph from popgraph.

[adj pop_labels] = pops_to_routers(popgraph,rules);