%test_sample_from_k4_dist

opts.vmin=10;
opts.vmax=100;
v_out=sample_from_k4_dist(1e7,opts);

vmag=vecnorm(v_out,2,2);

%%
sout=smooth_hist(vmag)
plot(sout.bin.centers,(sout.bin.centers.^4).*sout.count_rate.smooth)


%%
sout=smooth_hist(vmag,'scale_y_fun',@(x,y) y.*(x.^4),'sigma',1)
% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
plot(sout.bin.centers,sout.count_rate.smooth)
xlim([opts.vmin,opts.vmax])
%xlim([opts.vmin,opts.vmax])