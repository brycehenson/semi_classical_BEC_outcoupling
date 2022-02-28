%% compare the density profiles with experiment

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path
hebec_constants
rng('shuffle')

%%

%define exp parameters
% we will usea a 10mm/s oscillation because thats about what the tune out used
tf_param=[];
tf_param.trap_freq=2*pi*[51,412,415]; 
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
tf_param.num=6e5;
tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);

sim_opt=[];
sim_opt.osc.omega =tf_details.inputs.omega;
%lets define the amplitude in velocity then convert to spatial
osc_amp_vel=col_vec([0,10,0.00]*1e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.lambda=0.30;%5/500;
sim_opt.osc.phase=col_vec([1/4,0,(0)])*pi; %col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.grav=const.g0;
sim_opt.nsamp=8e5;%number of atoms to simulate
sim_opt.osc.time_offset=0;
sim_opt.t_end=12e-3;
sim_opt.verbose=2;


% do outoupling simulation
sim_out=class_sim_outcoupling(tf_details,sim_opt);


%%


fh=stfig('den hist');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
    
% plots
hist_opt=[];
hist_opt.cords_to_plot=[2,3];
hist_opt.cords_labels={'x','y','z'};
hist_opt.cords_units=repmat({'m/s'},[1,3]);
hist_opt.mean_marker=true;
hist_opt.disp_scale=[1,1,1]*1e-3;

hist_opt.kernel_factor=10;
% opts.std_norm_bounds=true;
% opts.hist_bounds=[5,1,1]*3;
hist_opt.std_norm_bounds=false;
%opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
hist_opt.cmap=viridis(512);
hist_opt.density_power=0.4;
hist_opt.hist_bounds=[[-50,50];[-40,40];[-30,80]]*1.1*1e-3;
hist_opt.blur_size=[1,1,1]*2.5e-4; %0.3e-3
% optional
hist_opt.filter_in_other_axis=false;
hist_opt.norm_den_to_num_in=true;
hist_opt.scale_den_disp_fac=true;
hist_opt.norm_den_unity=true;
%opts.cbar_lims=[0,2e-3];
% opts.bin_factor;
hist_opt.save_hist_as_image=false;
hist_opt.open_fig=false;
%hist_opt.save_hist_name=strcat(save_path,sprintf('raw_%0.4u.png',jj));


hist_fig=hist_2d_data( sim_out.final.vel.vals,hist_opt);


hold on
plot(sim_out.bec.vel(hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
    sim_out.bec.vel(hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
    '^','MarkerSize',10,'LineWidth',1.4,'Color',[0.1,0.1,1])
% plot(sim_out.dyn_states.bec.vel(time_idxs(ii),hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
%     sim_out.dyn_states.bec.vel(time_idxs(ii),hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
%     'wo','MarkerSize',5,'LineWidth',1.3)
hold off




set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=500;
fig_aspect_ratio=1; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.04,0])
%grid on
fig_name='pal_den_plot';
fig_dir='./figs/thesis_figs';
% export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
% %export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
% exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))


%% 1d density profile in y


colors_main=[[210,56,46];[50,163,73];[73,133,173]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);

xscale=1e3;

fh=stfig('1d den ');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
cord_to_plot=2;
    

sim_den=smooth_hist(sim_out.final.vel.vals(:,cord_to_plot)-sim_out.final.vel.mean(cord_to_plot),...
    'lims',hist_opt.hist_bounds(1,:),'sigma',2.5e-4);
norm_den=sim_den.count_rate.smooth/max(sim_den.count_rate.smooth);
plot( sim_den.bin.centers*xscale,norm_den,'color','k','LineWidth',1.5)


ylim([-0.04,1.1])
area_under_curve_sim=sum(norm_den)*mean(diff(sim_den.bin.centers));

x_samp=linspace(min(sim_den.bin.centers),max(sim_den.bin.centers),1e2)';
y_top_hat=abs(x_samp)<tf_details.pal_vmax;
area_under_curve_tophat=sum(y_top_hat)*mean(diff(x_samp));
y_top_hat=y_top_hat*(area_under_curve_sim/area_under_curve_tophat);

hold on
plot(x_samp*xscale,y_top_hat,':','LineWidth',1.5,'color',colors_main(2,:))


% now lets make a sampling distibution

% now lets make a sampling distibution from the experimenal distribution
samp_repeats=1e6;
pop_size=10;
sim_dist_mean_reps=dist_of_sample_means_discreet(sim_den.count_rate.smooth,sim_den.bin.centers,pop_size,samp_repeats);

sigma_smooth=1e-4;
%*sqrt(pop_size)
sim_den_mean=smooth_hist(sim_dist_mean_reps,'lims',hist_opt.hist_bounds(1,:),'sigma',sigma_smooth);
area_under_curve_sim=area_under_curve_sim*0.32;
area_under_curve=sum(sim_den_mean.count_rate.smooth)*mean(diff(sim_den_mean.bin.centers));
norm_den=sim_den_mean.count_rate.smooth*area_under_curve_sim/area_under_curve;
plot( sim_den_mean.bin.centers*xscale,norm_den,'--','color',colors_main(3,:),'LineWidth',1.5)
xlabel('$v_y-\bar{v_y}$ (mm/s)')
ylabel('density (mm/s)')
ylim([-0.04,1.1])


% now lets make a sampling distibution from the uniform dist
top_hat_samp_mean_reps=dist_of_sample_means_discreet(y_top_hat,x_samp,pop_size,samp_repeats);
%*sqrt(pop_size)
mean_top_sh=smooth_hist(top_hat_samp_mean_reps,'lims',hist_opt.hist_bounds(1,:),'sigma',sigma_smooth);
area_under_curve=sum(mean_top_sh.count_rate.smooth)*mean(diff(mean_top_sh.bin.centers));
norm_den=mean_top_sh.count_rate.smooth*area_under_curve_sim/area_under_curve;
plot( mean_top_sh.bin.centers*xscale,norm_den,'-.','color',colors_main(1,:),'LineWidth',1.5)
xlabel('mean subtracted velocity,$v_y-\bar{v_y}$ (mm/s)')
ylabel('density (arb. u.)')
ylim([-0.04,1.1])
hold off

ln=legend('simulation','uniform approx.','sample mean sim.','sample mean uniform')
legend('boxoff')
legend('location','northwest')
ln.EdgeColor=[1,1,1];

set(gca, {'XColor', 'YColor'}, {'k','k'});
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=500;
fig_aspect_ratio=0.8; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])
set(gca,'linewidth', 1.3)
set(gca,'TickLength',[0.02,0])
%grid on
fig_name='pal_1d_den_plot';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))

%%


