function pts=sample_pts_from_tf(n_sample,tf_details)
% sample some points from a 3d tf distribution using rejection sampling
% there are better ways to do sampling from a distribution but this is the easiest to implement

% set the number of pts to generate in each loop of the rejection sampling
% limited by memory
max_chunk_size=3e6;
verbose=0;

expected_sucess_rate=0.208; % this is independent of the trap freqs


pts=nan(n_sample,3);
pts_done=0;
pts_remaining=n_sample;

pad_factor=1.001; % pad some variables to make sure the edges are covered
cut_extra=true; %remove extra points from the bootstrap so that output is excatly n_sample long

while pts_remaining>0
    % set a chunk size, this will give the right number of sucess. sample pts on average
    % there may be some opt that could be done here for instance adding sqrt(chunk_size) to chunk_size
    % so that not sampling enough is 1sd down in prob
    chunk_size=ceil(pts_remaining/expected_sucess_rate);
    chunk_size=min([chunk_size,max_chunk_size]);
    % to start with we generate some points in a 3d unit box
    pos_sample=rand(chunk_size,3);
    pos_sample=2*(pos_sample-0.5); % make into random number on the -1,1 interval
    % now scale by the thomas fermi radius times a small padd factor
    pos_sample=pos_sample.*repmat(tf_details.tf_radi',[chunk_size,1])*pad_factor;
    % find what the tf density is at these locations
    density=tf_density_distribution(pos_sample,tf_details);
    % to do rejection sampling we generate another random number between 0 and the peak density
    amplitude_rand=rand(chunk_size,1)*tf_details.density_peak*pad_factor;
    % the sampled point is kept if the sample is less than the density
    mask=amplitude_rand<density;
    pts_tmp=pos_sample(mask,:);
    npts_chunk=size(pts_tmp,1);
    if npts_chunk>0
        pts(pts_done+1:pts_done+npts_chunk,:)=pts_tmp;
    end
    pts_done=pts_done+npts_chunk;
    pts_remaining=pts_remaining-npts_chunk;
    if verbose>1
        fprintf('chunk size %.2e sample sucess rate %f \n',chunk_size,npts_chunk/chunk_size)
    end
end


if cut_extra
    pts=pts(1:n_sample,:);
end

end