function sample_means=dist_of_sample_means_discreet(prob,cord,pop_size,repeats)
% make a number of samples of the mean of pop_size drawn from the
% prob with cordinate cord
sample_means=nan(repeats,1);
size_sim_den_norm=numel(prob);
parfor ii=1:repeats
    idxs_tmp=randsample(size_sim_den_norm,pop_size,true,prob);
    sample_means(ii)=mean(cord(idxs_tmp));
end

end