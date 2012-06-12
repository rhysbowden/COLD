function indices = random_subset(set_size,subset_size)
% indices = random_subset(set_size,subset_size)
% generates a uniform random subset of size 'subset_size' of the first
% 'set_size' integers. output is a row vector.
if(subset_size>set_size)
    indices = [];
    fprintf(2,'subset size can''t be greater than set size in random_subset');
else
    if(subset_size>0.5*set_size)
        actual_size=set_size-subset_size;
    else
        actual_size = subset_size;
    end
    
    index_array = zeros(1,set_size);
    for i=0:(actual_size-1)
        unchosen = find(index_array==0);
        index_array(unchosen(ceil(rand*(set_size-i))))=1;
    end

    if(subset_size>0.5*set_size)
        index_array = 1-index_array;
    end
    indices = find(index_array);
end
end