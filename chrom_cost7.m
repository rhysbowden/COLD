% chrom_cost7.m
% calculates the cost of a chromosome in top_gen.
% chrom_cost2: updated version of chrom cost to include k4
% chrom_cost3: changes the input and output formats.
% chrom_cost4: just the overall city link distances are supplied.
% chrom_cost5: works just like chrom_cost4, but doesn't break for Inf
% distances. (About 10% slower than chrom_cost4 ??)
% chrom_cost7: changed to work with disconnected adjacency matrices -
% different input/output now too. Use chrom_cost5 or 6 for adj matrices
% guaranteed to be connected.
function [total_cost chromosome_changed new_chromosome cost_breakdown] = chrom_cost7(A,distances,demands,kparams,distance_order)
A = A>0;
new_chromosome = A;
chromosome_changed=0;
n = size(A,1);
k0 = kparams(1);
k1 = kparams(2);
k2 = kparams(3);
k4 = kparams(4);
BIG = n*(max(max(distances))+1);
% chrom dists is the normal city distances for links in the chromosome, and
% a large number times the rank order of the link if it's not in the
% chromosome.
chrom_dists = zeros(n);
chrom_dists(A>0) = distances(A>0);
chrom_dists(A==0) = BIG*distance_order(A==0);
% do the routing
[path_distances pred] = all_shortest_paths(sparse(chrom_dists));
if(any(any(path_distances>=BIG)))
    chromosome_changed = 1;
    new_chromosome = zeros(n);
    % find the links that are used.
    for column =1:n
        new_chromosome(pred(pred(:,column)>0,column),column)=1;
    end
    A = new_chromosome;
    chrom_dists = zeros(n);
    chrom_dists(A>0) = distances(A>0);
    chrom_dists(A==0) = BIG*distance_order(A==0);
    [path_distances] = all_shortest_paths(sparse(chrom_dists));
end

% find the costs
k2cost = k2*sum(sum(path_distances(demands>0).*demands(demands>0)))/2;
k1cost = k1*sum(sum(distances(A>0.5)))/2;
k0cost = k0*sum(sum(A))/2;
k4cost = k4*sum(sum(A)>1);
total_cost = k2cost+k1cost+k0cost+k4cost;
cost_breakdown = [k0cost k1cost k2cost k4cost];
end