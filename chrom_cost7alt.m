% chrom_cost7alt.m
% calculates the cost of a chromosome in top_gen.
% chrom_cost2: updated version of chrom cost to include k3
% chrom_cost3: changes the input and output formats.
% chrom_cost4: just the overall city link distances are supplied.
% chrom_cost5: works just like chrom_cost4, but doesn't break for Inf
% distances. (About 10% slower than chrom_cost4 ??)
% chrom_cost7: changed to work with disconnected adjacency matrices -
% different input/output now too. Use chrom_cost5 or 6 for adj matrices
% guaranteed to be connected.

% modified with Yin Zhang's Floyd-Warshall algorithm script, rather than
% MatlabBGL's; a little faster but requires MEX compilation
function [total_cost,chromosome_changed,new_chromosome,cost_breakdown] = chrom_cost7alt(A,distances,demands,kparams,distance_order)
A = A>0;
new_chromosome = A;
chromosome_changed=0;
n = size(A,1);
k0 = kparams(1);
k1 = kparams(2);
k2 = kparams(3);
k3 = kparams(4);
BIG = n*(max(max(distances))+1);
% chrom dists is the normal city distances for links in the chromosome, and
% a large number times the rank order of the link if it's not in the
% chromosome.
chrom_dists = zeros(n);
chrom_dists(A>0) = distances(A>0);
chrom_dists(A==0) = BIG*distance_order(A==0);
% do the routing
[path_distances, pred] = floyd_apsp(chrom_dists);
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
    [path_distances] = floyd_apsp(chrom_dists);
end



[~,~,K] = size(demands);
% this simple method is fine, all costs are non-negative
k2cost = 0; % minimax criterion; defaults to single TM if only one TM in use
for i=1:K
    demand = demands(:,:,i);
    k2cost = max(k2cost,k2*sum(sum(path_distances(demand>0).*demand(demand>0)))/2);
end

% find the costs
k0cost = k0*sum(sum(A))/2;
k1cost = k1*sum(sum(distances(A>0.5)))/2;
k3cost = k3*sum(sum(A)>1);
total_cost = k2cost+k1cost+k0cost+k3cost;
cost_breakdown = [k0cost k1cost k2cost k3cost];
end