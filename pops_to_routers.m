function [router_graph pop_labels] = pops_to_routers(popgraph,rules)
% [router_graph pop_labels] = pops_to_routers(popgraph,rules)
% takes a PoP level network topology with n PoPs and produces a router level
% network topology using generalised graph products.
% 
% input:
% popgraph is a struct with: 
% popgraph.adjacency = the n by n adjacency matrix
% popgraph.traffic = an n by n matrix with popgraph.traffic(i,j) = the total
%   traffic along link (i,j), 0 if (i,j) doesn't exist.
% popgraph.access_traffic = an n by 1 vector of access traffics for each
%   node, that is, traffic that does not transit the PoP but terminates
%   there.
% popgraph.node_types = vector denoting the type of each node, used for type
%   rules.
% popgraph.num_access_routers = the number of access routers required in
% each PoP.
%
% rules is a struct with:
% rules.templates = a row cell array with templates{i} the H fiber for node
%   i; that is, the core router topology at PoP i.
% rules.products = a 2 dimensional cell array with rules.templates{i,j} 
%   = the template for 
%   the links from type i nodes to type j nodes. each template is an 
%   adjacency matrix. for cartesian products use eye(sizei,sizej) and for
%   strong product use ones(sizei,sizej) where sizei is the size of template
%   i.
% rules.pop_discriminator = 'traffic', 'type' or 'hubs'; determines the way router graphs
%   are allocated to PoPs.
% rules.link_discriminator = 'traffic' or 'ends'; determines the way links
%   are determined between PoPs, either by the PoPs at the end or how much
%   traffic is flowing down the link.
% rules.thresholds = a vector of the traffic thresholds. 
%   If rules.discriminator = 'traffic' then all nodes with
%   traffic less than thresholds(1) will be type 1, all nodes with traffic
%   between threshold(1) and threshold(2) will be type 2, etc.
% rules.access_thresholds = an n by 1 vector of thresholds for access 
%   traffic that determines how many access routers are required at each PoP.
% rules.link_thresholds = a threshold for changing between cartesian
%   product and strong product
%
% output:
% router_graph = adjacency matrix of the router level graph.
% pop_labels = a vector listing the number of the pop to which each router 
%   belongs. pop_labels(i) = the pop containing router i.
%
% see pops_to_routers_example.m for a simple example of using this function
% with default parameters.

% useful variables
num_pops = size(popgraph.adjacency,1);

% compulsory fields
adjacency = popgraph.adjacency;
% defaults
pop_discriminator = 'hubs';
link_discriminator = 'ends';
node_types = (sum(adjacency)>1)+1;
templates = {0,[0 1; 1 0]};
products = {1,[1 1];[1;1],[1 0; 0 1]};
%num_access_routers = 2*ones(1,num_pops);
% optional fields
if(isfield(popgraph,'traffic'))
    traffic = popgraph.traffic;
end
if(isfield(popgraph,'node_types'))
    node_types = popgraph.node_types;
end
if(isfield(rules,'templates'))
    products = rules.templates;
end
if(isfield(rules,'products'))
    products = rules.products;
end
if(isfield(rules,'pop_discriminator'))
    pop_discriminator = rules.pop_discriminator;
end
if(isfield(rules,'link_discriminator'))
    link_discriminator = rules.link_discriminator;
end
if(isfield(rules,'thresholds'))
    thresholds = rules.thresholds;
end
if(isfield(rules,'access_thresholds'))
    access_thresholds = rules.access_thresholds;
end
if(isfield(popgraph,'num_access_routers'))
    num_access_routers = popgraph.num_access_routers;
    if(length(popgraph.num_access_routers)~=size(popgraph.adjacency,1))
        throw(MException('COLD:InvalidInput','popgraph.num_access_routers and popgraph.adjacency have differing numbers of routers.'));
    end
else
    if(isfield(rules,'access_thresholds'))
        access_thresholds = rules.access_thresholds;
    else
        throw(MException('COLD:InvalidInput','rules.access_thresholds must be defined if popgraph.num_access_routers is not.'));
    end
    if(isfield(popgraph,'access_traffic'))
        access_traffic = popgraph.access_traffic;
    else
        throw(MException('COLD:InvalidInput','rules.access_traffic must be defined if popgraph.num_access_routers is not.'));
    end
end

% find the node types based on their traffic and the traffic thresholds.
if(strcmp(pop_discriminator,'hubs'))
    node_types = ones(1,num_pops)+ (sum(adjacency)>1);
elseif(strcmp(pop_discriminator,'traffic'))
    in_traffic = sum(traffic,1);
    out_traffic = sum(traffic,2)';
    node_traffic = in_traffic+out_traffic; % transit traffic is counted twice, terminal traffic is only counted once.
    node_types = max((repmat(node_traffic,length(thresholds)+1,1)>repmat([-Inf thresholds]',1,length(node_traffic))).*repmat((1:(length(thresholds)+1))',1,length(node_traffic)));
elseif(strcmp(pop_discriminator,'types'))
else
    throw(MException('COLD:InvalidInput','popgraph.pop_discriminator should be ''hubs'',''traffic'' or ''types'''));
end
% find the number of access routers for each node
if(~isfield(popgraph,'num_access_routers'))
    num_access_routers = max((repmat(access_traffic,length(access_thresholds)+1,1)>repmat([-Inf access_thresholds]',1,length(access_traffic))).*repmat((1:(length(access_thresholds)+1))',1,length(access_traffic)));
end
% get the size of each template so we can start building the output
% router_graph
template_sizes = zeros(1,size(products,2));
for i = 1:size(products,2)
    template_sizes(i) = size(products{i,i},2);
end
node_sizes = template_sizes(node_types);
router_graph = zeros(sum(node_sizes)+sum(num_access_routers),sum(node_sizes)+sum(num_access_routers)); % the output
pop_labels = zeros(1,size(router_graph,1));
start_indices = (cumsum(node_sizes)-node_sizes)+1; % the first row in the matrix corresponding to node i.
end_indices = cumsum(node_sizes);

% form the core network
for i = 1:num_pops
    pop_labels(start_indices(i):end_indices(i))=i;
    for j = 1:num_pops
        if(i==j)
            router_graph(start_indices(i):end_indices(i),start_indices(i):end_indices(i)) = templates{node_types(i)};
        else
            if(strcmp(link_discriminator,'ends'))
                router_graph(start_indices(i):end_indices(i),start_indices(j):end_indices(j)) = adjacency(i,j)*products{node_types(i),node_types(j)};
            elseif(strcmp(link_discriminator,'traffic'))
                % cartesian product used if traffic low enough, otherwise
                % strong product
                if(traffic(i,j)<=link_thresholds(1))% cartesian product
                    router_graph(start_indices(i):end_indices(i),start_indices(j):end_indices(j)) = adjacency(i,j)*eye(node_sizes(i),node_sizes(j));
                else % strong product
                    router_graph(start_indices(i):end_indices(i),start_indices(j):end_indices(j)) = adjacency(i,j)*ones(node_sizes(i),node_sizes(j));
                end
            else
                throw(MException('COLD:InvalidInput','link_discriminator should be ''traffic'' or ''ends'''));
            end
        end
    end
end
% form the access router network
start_indicesAR = (cumsum(num_access_routers)-num_access_routers)+1+sum(node_sizes); % the first row in the matrix corresponding to node i.
end_indicesAR = cumsum(num_access_routers)+sum(node_sizes);
for i=1:num_pops
    pop_labels(start_indicesAR(i):end_indicesAR(i))=i;
    router_graph(start_indices(i):end_indices(i),start_indicesAR(i):end_indicesAR(i)) = ones(node_sizes(i),num_access_routers(i));
    router_graph(start_indicesAR(i):end_indicesAR(i),start_indices(i):end_indices(i)) = ones(num_access_routers(i),node_sizes(i));
end
% todo: add node_positions and ability to create a graph