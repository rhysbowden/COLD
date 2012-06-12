function topology = coldGA(parameters)
% coldGA:
% Runs a genetic algorithm to create a model of an Internet topology
% See below for explanation of parameters.
% Rhys Bowden Sep 2009
% version 11 - first github version
% Feb 2011
% requires matlab_bgl package for all_shortest_paths
% http://dgleich.github.com/matlab-bgl/

display(' --- topology generator 11 ---');
topology = struct();
topology.version = [11 0]; % version 11.0 of the generator

% ---- parameters ---- 
% choose positive parameters (except k1 which can be zero, and k2 or k0)
% number of nodes
num_nodes = parameters.num_nodes;
% width of zone
width = parameters.width;
% height of zone
height = parameters.height;
% population parameter (can change population distribution later)
node_population_param = parameters.node_population_param;
% number of chromosomes in population
num_chromosomes=parameters.num_chromosomes;
% number of generations in the genetic algorithm
num_generations=parameters.num_generations;
% optimisation constants
k2 = parameters.k2;% bandwidth*length constant
k1 = parameters.k1;% length constant
k0 = parameters.k0;% existence constant
k4 = parameters.k4; % multi-connection penalty %NEW to version 8
% crossover parameters
% we assign these to a and b later, in the crossover section
crossovera=parameters.crossovera;
crossoverb=parameters.crossoverb;
% probability of mutation
mutate_fn = parameters.mutate_fn; % function handle. called here with no parameters.
% average number of links per node.
start_links = parameters.start_links;
% number of chromosomes to save from each generation
num_saved_chromosomes = min(parameters.num_saved_chromosomes,parameters.num_chromosomes);
if(ismember('pop_distn',fieldnames(parameters)))
    pop_distn = parameters.pop_distn;
else
    pop_distn = 'exp';
end

% ----- random initial conditions ----
% generate the node populations and positions if not given as input
% use the positions and populations to get pairwise distances and traffic
% matrix.
if(ismember('node_distances',fieldnames(parameters)))
    node_distances = parameters.node_distances;
    demand = parameters.demand;
else
    node_map_determined_by_function = 0;
    if(ismember('node_map_function',fieldnames(parameters)))
        node_map_function = parameters.node_map_function;
        node_map_determined_by_function = 1;
    end
    % -------- node map -----------
    % this part creates a node map array with dimensions (numnodes by 3)
    % each row is a node, the first column is the x coordinate, the
    % second the y coordinate and the third the population.
    if(~node_map_determined_by_function)
        placesX = unifrnd(0,width,[num_nodes 1]);
        placesY = unifrnd(0,height,[num_nodes 1]);
    else
        places = node_map_function(num_nodes);
        placesX = places(1:num_nodes,1);
        placesY = places(1:num_nodes,2);
    end

    if(strcmp(pop_distn,'pareto'))
        % in this case, node_population_param(1) is scale (xm), and
        % node_population_param(2) is shape (alpha), see wikipedia.
        pop_xm = node_population_param(1);
        pop_alpha = node_population_param(2);
        places_pop = unifrnd(0,1,[num_nodes 1])+eps;
        places_pop = pop_xm./(places_pop.^(1/pop_alpha));
    elseif(strcmp(pop_distn,'gamma'))
        places_pop=random(pop_distn,node_population_param(1),node_population_param(2),num_nodes,1);
    else
        places_pop = random(pop_distn,node_population_param,num_nodes,1);
    end
    node_map = [placesX placesY places_pop];

    % this generates a matrix of distances between nodes
    % and a parameter proportional to traffic demand (gravity model)
    node_distances = zeros(num_nodes);
    demand=zeros(num_nodes);
    for i=1:num_nodes
        for j=1:i-1
            distance_temp = sqrt((node_map(i,1)-node_map(j,1))^2+(node_map(i,2)-node_map(j,2))^2);
            node_distances(i,j) = distance_temp;
            node_distances(j,i) = distance_temp;
            demand_temp = node_map(i,3)*node_map(j,3);
            demand(i,j)=demand_temp;
            demand(j,i)=demand_temp;
        end
    end
end % end of distance and demand setup

% each potential link is assigned a number from 1 to n^2, in the same order
% as the link lengths. This is for implicit minimum spanning tree for
% disconnected networks.
[distance_order ix] = sort(reshape(node_distances,[],1));
distance_order=zeros(size(ix));
for ind = 1:length(ix)
    distance_order(ix(ind))=ind;
end
distance_order = reshape(distance_order,num_nodes,num_nodes);

% ----------- generate initial chromosomes ------------
% each chromosome is a matrix A with A(i,j)=1 iff there is a link from node i to node j
% chromosomes are adjacency matrices, stored in a 3d array, where each chromosome is a slice (:,:,i).
% as of this moment, initial chromosomes are generated randomly and
% uniformly
new_chromosomes = zeros(num_nodes,num_nodes,num_chromosomes);
% seed with minimum spanning tree and fully connected graph.
new_chromosomes(:,:,1) = full(mst(sparse(node_distances)));
if(num_chromosomes>=2)
new_chromosomes(:,:,2) = ones(size(node_distances))-eye(size(node_distances));
end
% if there are starting chromosomes given in the parameters use as many as
% possible
num_init = 0;
if(ismember('initial_chromosomes',fieldnames(parameters)))
    initial_chromosomes = parameters.initial_chromosomes;
    num_init = size(initial_chromosomes,3);
    fprintf(1,'NUMBER INITIAL CHROMOSOMES %d\n',num_init);
    new_chromosomes(:,:,3:min(num_chromosomes,num_init+2)) = initial_chromosomes(:,:,1:min(num_chromosomes-2,end));
end
% fill the rest with random chromosomes
new_chromosomes(:,:,(min(num_chromosomes,num_init+2)+1):num_chromosomes) = rand(num_nodes,num_nodes,num_chromosomes-min(num_chromosomes,num_init+2))<(start_links/num_nodes);
% make them symmetrical with zeros on the diagonal
for i=1:num_nodes
    for j=1:i-1
        new_chromosomes(j,i,:)=new_chromosomes(i,j,:);
    end
    new_chromosomes(i,i,:)=0;
end

% ---------- find the initial costs --------------------
new_chromosome_costs = zeros(1,num_chromosomes);
for i=1:num_chromosomes
    [new_chromosome_costs(i) chromosome_changed replacement] = chrom_cost7(new_chromosomes(:,:,i),node_distances,demand,[k0 k1 k2 k4],distance_order);
    if(chromosome_changed)
        % if the chromosome is disconnected this will be the minimum
        % (link-length-wise) connected version of it
        new_chromosomes(:,:,i) = replacement;
    end
end

% ----------- run the GA -------------------------------
for generation = 1:num_generations
    % --------- initialise -----------
    % chromosomes are the chromosomes from last generation.
    chromosomes = new_chromosomes;
    chromosome_costs = new_chromosome_costs;
    % new chromosomes for the next generation
    new_chromosomes = zeros(num_nodes,num_nodes,num_chromosomes);
    new_chromosome_costs = zeros(1,num_chromosomes);

    % --------- print costs -----------
    fprintf('Generation: %2d  Min cost: %7g Max cost: %7g Mean cost: %7g Cost SD: %g \n',generation,min(chromosome_costs),max(chromosome_costs),mean(chromosome_costs),sqrt(var(chromosome_costs))); % display!

    % ---------- mutation and crossover ------------------
    [temp ix] = sort(chromosome_costs,'ascend');
    new_chromosomes(:,:,1:num_saved_chromosomes) = chromosomes(:,:,ix(1:num_saved_chromosomes));
    new_chromosome_costs(1:num_saved_chromosomes) = chromosome_costs(ix(1:num_saved_chromosomes));
    % pick b at random, then choose the best a of them.
    for i=(num_saved_chromosomes+1):num_chromosomes
        if(rand<parameters.crossover_rate)
            % --- CROSSOVER ----
            indices = (random_subset(num_chromosomes,crossoverb))'; % the b chosen indices
            subcosts = chromosome_costs(indices);
            % pick best a out of b random (for the crossover)
            [temp sorted_indices] = sort(subcosts,'ascend');
            chosen_indices = indices(sorted_indices(1:crossovera)); % list of indices.

            % chosen chromosomes for crossover (a of them)
            chosen_costs = chromosome_costs(chosen_indices);%a*1
            chosen_costs_inverse=1./chosen_costs; %costs should be positive
            chosen_costs_inverse = chosen_costs_inverse/sum(chosen_costs_inverse); %normalise
            cum_prob = zeros(num_nodes);
            for cc_ind = 1:length(chosen_indices)
                cum_prob = cum_prob + chromosomes(:,:,chosen_indices(cc_ind))*chosen_costs_inverse(cc_ind);
            end
            new_chromosome = triu(rand(num_nodes)<cum_prob);
            new_chromosome = new_chromosome+new_chromosome';        

        else % mutate. mutation and crossover are complementary events now.
            cum_sum_costs = cumsum(chromosome_costs);
            new_chromosome = chromosomes(:,:,find(cum_sum_costs>rand*cum_sum_costs(end),1,'first'));
            % --- MUTATION ---
            hubs = find(sum(new_chromosome)>1);
            if(rand<parameters.hub_mut_rate&&length(hubs)>1) % hub mutation
                hubs = find(sum(new_chromosome)>1);

                hub_removed = hubs(ceil(rand*length(hubs)));
                new_hubs = sum(new_chromosome)>1;
                new_hubs(hub_removed)=0;
                new_hubs = find(new_hubs);
                [ctemp hub_link] = min(node_distances(new_hubs,hub_removed));
                new_row = zeros(1,num_nodes);
                new_row(new_hubs(hub_link))=1;
                new_chromosome(hub_removed,:) = new_row;
                new_chromosome(:,hub_removed) = new_row';
                new_disconnected = find(sum(new_chromosome)==0);
                if(length(new_hubs)>1)
                    [ctemp new_attachments] = min(node_distances(new_hubs,new_disconnected));
                else
                    new_attachments = new_hubs*ones(1,length(new_disconnected));
                end
                for attach_index = 1:length(new_disconnected)
                    new_chromosome(new_attachments(attach_index),new_disconnected(attach_index))=1;
                    new_chromosome(new_disconnected(attach_index),new_attachments(attach_index))=1;
                end
            else % link by link mutation
                mutations = mutate_fn();
                existent_links = find(triu(new_chromosome,1));
                nonexistent_links = find(triu(1-new_chromosome,1));
                num_existent = length(existent_links);
                num_nonexistent = length(nonexistent_links);
                mutations = min(mutations,[num_existent num_nonexistent]);
                neg_mut_ind = random_subset(num_existent,mutations(1));
                pos_mut_ind = random_subset(num_nonexistent,mutations(2));
                neg_mut = existent_links(neg_mut_ind);
                pos_mut = nonexistent_links(pos_mut_ind);
                new_chromosome(neg_mut) = 0;
                new_chromosome(pos_mut) = 1;
                % fix it so that it's symmetric with 0s on the diagonal.
                new_chromosome = triu(new_chromosome,1);
                new_chromosome = new_chromosome+new_chromosome';
            end
        end
        %add to the new population
        new_chromosomes(:,:,i)=new_chromosome;
    end % of all crossovers
    % ------- find the cost values for the new chromosomes -------
    for i=(num_saved_chromosomes+1):num_chromosomes
        [new_chromosome_costs(i) chromosome_changed replacement] = chrom_cost7(new_chromosomes(:,:,i),node_distances,demand,[k0 k1 k2 k4],distance_order);
        if(chromosome_changed)
            % if the chromosome is disconnected this will be the minimum
            % (link-length-wise) connected version of it
            new_chromosomes(:,:,i) = replacement;
        end
    end
end % of all generations
[best_cost best_index]=min(new_chromosome_costs);
best_chromosome = new_chromosomes(:,:,best_index);
display(best_cost);

% ----- OUTPUTS -----
if(~ismember('node_distances',fieldnames(parameters)))
    topology.node_map = node_map;
end
topology.adjacency = best_chromosome;
topology.node_distances = node_distances;
topology.demand = demand;
% also outputs topolgy.version. See top.
end