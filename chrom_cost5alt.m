% chrom_cost5alt.m
% calculates the cost of a chromosome in top_gen.
% chrom_cost2: updated version of chrom cost to include k4
% chrom_cost3: changes the input and output formats.
% chrom_cost4: just the overall city link distances are supplied.
% chrom_cost5: works just like chrom_cost4, but doesn't break for Inf
% distances. (About 10% slower than chrom_cost4 ??)

% modified with Yin Zhang's Floyd-Warshall algorithm script, rather than
% MatlabBGL's
function [total_cost,cost_breakdown,chrom_dists] = chrom_cost5alt(A,distances,demands,kparams)
A = A>0;
BIG = size(A,1)*(max(max(distances))+1);
% do the routing
chrom_dists = zeros(size(A));
chrom_dists(A>0) = distances(A>0);
chrom_dists(A==0) = BIG;
path_distances = floyd_apsp(chrom_dists);
% find the costs
k4cost = kparams(4)*sum(sum(A)>1);
k2cost = kparams(3)*sum(sum(path_distances(demands>0).*demands(demands>0)))/2;
k1cost = kparams(2)*sum(sum(distances(A>0.5)))/2;
k0cost = kparams(1)*sum(sum(A))/2;
total_cost = k2cost+k1cost+k0cost+k4cost;
cost_breakdown = [k0cost k1cost k2cost k4cost];
end