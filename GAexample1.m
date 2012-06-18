% an example of setting parameters and calling the Genetic Algorithm part
% of COLD.
% see github.com/rhysbowden/COLD
num_nodes = 40;
width =10;
height = 10;
node_population_param = 30;
pop_distn = 'exp';
% optimisation parameters
parameters = struct();
parameters.k0 = 10;
parameters.k1 = 1;
parameters.k2 = 3e-4;
parameters.k4 = 50;
% internal GA parameters
parameters.num_nodes = num_nodes;
parameters.width = width;
parameters.height = height;
parameters.node_population_param = node_population_param;
parameters.pop_distn = 'exp';
parameters.num_chromosomes=100;
parameters.num_generations=100;
parameters.crossovera=2;
parameters.crossoverb=10;
parameters.mutate_prob = 0.01;
parameters.start_links = 2;
parameters.mutate_fn = @()geornd([0.4,0.4]);
parameters.num_saved_chromosomes = 50;
parameters.fast_mode =1;
parameters.hub_mut_rate = 0.05;
parameters.crossover_rate = 0.3;

% the setup for the optimisation
placesX = unifrnd(0,width,[num_nodes 1]);
placesY = unifrnd(0,height,[num_nodes 1]);
places_pop = random(pop_distn,node_population_param,num_nodes,1);
node_map = [placesX placesY places_pop];
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
parameters.demand = demand;
parameters.node_distances = node_distances;

topo = coldGA(parameters);