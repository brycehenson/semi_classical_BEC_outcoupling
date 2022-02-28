function v_out=sample_from_k4_dist(n_sample,opts)


vmin=opts.vmin;
vmax=opts.vmax;

rand_dir= randn(n_sample,3); % random direction
% normalize each to be on the unit sphere
rand_dir=rand_dir./repmat(vecnorm(rand_dir,2,2),1,3); 

% to sample the v^(-4) distribution f(v)=c v^-4
% we sample in the transformed cordinate q=v^-3
% v=q^(-1/3)
% which has the prob dist f_q (q)= c/3
% which is fantastic and uniform
% and then transform the sampled value back to the velocity space


q_samp=rand(n_sample,1);
% now we must scale the sample to the limits porvided
qmax=vmin^(-3);
qmin=vmax^(-3);
q_samp=qmin+q_samp*(qmax-qmin);
% now we transform to the velocity
v_samp=q_samp.^(-1/3);

v_out=rand_dir.*repmat(v_samp,1,3);


end