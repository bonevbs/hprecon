% generate cluster tree using bisection
function cluster = gen_cluster_rec(cluster, block_size)
  if cluster(1) > block_size || cluster(2) - cluster(1) > block_size
    if cluster(1) < block_size
      % refine only the second cluster
      cl2 = [floor((cluster(2)-cluster(1))/2), cluster(2) - cluster(1)];
      cl2 = gen_cluster_rec(cl2, block_size) + cluster(1);
      cl1 = cluster(1)*ones(1, length(cl2));
    elseif  cluster(2) - cluster(1) < block_size
      % refine only the first cluster
      cl1 = [floor(cluster(1)/2), cluster(1)];
      cl1 = gen_cluster_rec(cl1, block_size);
      cl2 = cluster(2)*ones(1, length(cl1));
    else
      % refine both
      cl1 = [floor(cluster(1)/2), cluster(1)];
      cl1 = gen_cluster_rec(cl1, block_size);
      cl2 = [floor((cluster(2)-cluster(1))/2), cluster(2) - cluster(1)];
      cl2 = gen_cluster_rec(cl2, block_size) + cluster(1);
      % make them equally long
      m = max(length(cl1), length(cl2));
      cl1 = ones(m/length(cl1),1) * cl1;
      cl2 = ones(m/length(cl2),1) * cl2;
    end
    cluster = [cl1(:)', cl2(:)'];
  end
end