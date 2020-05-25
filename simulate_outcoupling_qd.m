%simulate_outcoupling

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path

hebec_constants
rng('shuffle')
return

%%

%define exp parameters
param=[];
param.trap_freq=2*pi*[50,500,500];
param.num=1e6;
param.mass=const.mhe;
param.a_scat_len=const.ahe_scat;

tf_details=bec_properties(param.trap_freq/(2*pi),param.num,param.mass,param.a_scat_len)
%% use rejection sampling to sample from a 3d TF
n_sample=3e6;
tic 
den_power=2;
pos_sample=sample_pts_from_tf(n_sample,tf_details,den_power);
toc
mom_sample=1

%%
slice_size=0.1;
selection_mask=abs(pos_sample(:,2))<slice_size*tf_details.tf_radi(2) & abs(pos_sample(:,3))<slice_size*tf_details.tf_radi(3) ;
out_struct=smooth_hist(pos_sample(selection_mask,1),'lims',[-1,1]*tf_details.tf_radi(1),'sigma',tf_details.tf_radi(1)*1e-2,'bin_factor',10);
stfig('count_rate','add_stack',1); %this time we will give the figure a name and prepend the function that called it
clf
plot(out_struct.bin.centers,out_struct.count_rate.smooth*1e-3)
ylabel('count rate ($\mathrm{m}^{-1}$)')
xlabel('x pos (m)')

%%
opts=[];
opts.cords_to_plot=[2,3];
opts.cords_labels={'x','y','z'};
opts.cords_units={'m','m','m'};
opts.disp_scale=[1,1,1]*1e-6;
opts.blur_size=[1,1,1]*0.02;
opts.std_norm_bounds=true;
opts.hist_bounds=[5,1,1]*3;
% optional
opts.filter_in_other_axis=true;
% opts.bin_factor;

hist_2d_data(pos_sample,opts)

%% do a calssical sim of the QD atoms rolling of the BEC TF potential
sim_opt=[];
% sample from TF density squared
sim_opt.den_power=2;
% lets use a stationary BEC
sim_opt.osc.omega =tf_details.inputs.omega;
%lets define the amplitude in velocity then convert to spatial
sim_opt.osc.amp=col_vec([0,0,0]); 
sim_opt.osc.phase=col_vec([0,0,0]);
sim_opt.nsamp=3e5;
%sim_opt.grav=const.g0;
sim_opt.grav=0;
sim_opt.qd=[];
% bec is properly finished by 50e-3 a larger inital velociy allows a smaller t_end
%sim_opt.qd.vmin=50e-3;
sim_opt.qd.vmin=40e-3;
sim_opt.qd.vmax=400e-3;
sim_opt.t_end=0.03;

sim_out_bec_static=class_sim_outcoupling(tf_details,sim_opt);
%

%% lets look at the spatial distribution
opts=[];
opts.cords_to_plot=[2,3];
opts.cords_labels={'x','y','z'};
opts.cords_units=repmat({'m $\mathrm{s}^{-1}$'},[1,3]);
opts.mean_marker=true;
opts.disp_scale=[1,1,1]*1e-6;
opts.blur_size=[1,1,1]*0.02;
% opts.std_norm_bounds=true;
% opts.hist_bounds=[5,1,1]*3;
opts.std_norm_bounds=false;
opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
opts.hist_bounds=[1,1,1]*3;
opts.std_norm_bounds=true;
% optional
opts.filter_in_other_axis=false;
opts.norm_den_to_num_in=true;
opts.scale_den_disp_fac=true;
opts.norm_den_unity=true;
opts.save_hist_as_image=true;
%opts.save_hist_name='initial_pos_dist.tiff';

%opts.cbar_lims=[0,2e-3];
% opts.bin_factor;

%hist_2d_data(sim_out.start.vel.vals,opts)
hist_2d_data(sim_out_bec_static.start.pos.vals,opts)

hold on
plot(sim_out_bec_static.bec.vel(opts.cords_to_plot(1))/opts.disp_scale(opts.cords_to_plot(1)),...
    sim_out_bec_static.bec.vel(opts.cords_to_plot(2))/opts.disp_scale(opts.cords_to_plot(1)),...
    'wo','MarkerSize',10,'LineWidth',1.0)
hold off

% comapre the mean postion in velocity space to the bec
% this tells us how much the BEC gets pushed arround from the outcoupling process
% 70(7)e-6 ms⁻¹ is what i got in a 1e7 count sim
% the impulse on the bec is therefore 70e-6*nout*mass
% the change in velocity of the bec is therefor 70e-6*nout*mass=m*deltav*nin
% 70e-6 ms⁻¹*nout/nin=deltav
% 70e-6  ms⁻¹*(N_0*eta)/(N_0*(1-eta))=deltav
% for small outcoupling we aproximate this as 
% 70e-6  ms⁻¹*eta=deltav
% for a eta=1e-2
% \delta v=0.7 µs⁻¹


str_vals={};
disp_fac=1e6;
for ii=1:3
    str_vals{ii}=string_value_with_unc(sim_out_bec_static.final.vel.mean(ii)*disp_fac,...
        disp_fac*sim_out_bec_static.final.vel.std(ii)/sqrt(size(sim_out_bec_static.final.vel.vals,1)),'type','b');
end
fprintf('mean vel (x,y,z)= ( %s, %s, %s ) µms⁻¹ \n',str_vals{:})
%set(gcf,'Position',[100 100 1000 1000])
%%

slice_size=0.1;
pos_sample=sim_out_bec_static.start.pos.vals;
selection_mask=abs(pos_sample(:,2))<slice_size*tf_details.tf_radi(2) & abs(pos_sample(:,3))<slice_size*tf_details.tf_radi(3) ;
out_struct=smooth_hist(pos_sample(selection_mask,1),'lims',[-1,1]*tf_details.tf_radi(1),'sigma',tf_details.tf_radi(1)*1e-2,'bin_factor',10);
stfig('count_rate','add_stack',1); %this time we will give the figure a name and prepend the function that called it
clf
plot(out_struct.bin.centers,out_struct.count_rate.smooth*1e-3)
ylabel('count rate ($\mathrm{m}^{-1}$)')
xlabel('x pos (m)')

%% plots
opts=[];
opts.cords_to_plot=[2,3];
opts.cords_labels={'x','y','z'};
opts.cords_units=repmat({'m $\mathrm{s}^{-1}$'},[1,3]);
opts.mean_marker=true;
opts.disp_scale=[1,1,1]*1e-3;
opts.blur_size=[1,1,1]*1e-3;
% opts.std_norm_bounds=true;
% opts.hist_bounds=[5,1,1]*3;
opts.std_norm_bounds=false;
opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
opts.hist_bounds=[[-50,50];[-100,100];[-100,100]]*1e-3;
% optional
opts.filter_in_other_axis=false;
opts.norm_den_to_num_in=true;
opts.scale_den_disp_fac=true;
opts.norm_den_unity=true;
opts.save_hist_as_image=true;
%opts.save_hist_name='final_vel_dist.tiff';

%opts.cbar_lims=[0,2e-3];
% opts.bin_factor;

hist_2d_data(sim_out_bec_static.final.vel.vals,opts)

hold on
plot(sim_out_bec_static.bec.vel(opts.cords_to_plot(1))/opts.disp_scale(opts.cords_to_plot(1)),...
    sim_out_bec_static.bec.vel(opts.cords_to_plot(2))/opts.disp_scale(opts.cords_to_plot(1)),...
    'wo','MarkerSize',10,'LineWidth',1.0)
hold off

% comapre the mean postion in velocity space to the bec
% this tells us how much the BEC gets pushed arround from the outcoupling process
% 70(7)e-6 ms⁻¹ is what i got in a 1e7 count sim
% the impulse on the bec is therefore 70e-6*nout*mass
% the change in velocity of the bec is therefor 70e-6*nout*mass=m*deltav*nin
% 70e-6 ms⁻¹*nout/nin=deltav
% 70e-6  ms⁻¹*(N_0*eta)/(N_0*(1-eta))=deltav
% for small outcoupling we aproximate this as 
% 70e-6  ms⁻¹*eta=deltav
% for a eta=1e-2
% \delta v=0.7 µs⁻¹


str_vals={};
disp_fac=1e6;
for ii=1:3
    str_vals{ii}=string_value_with_unc(sim_out_bec_static.final.vel.mean(ii)*disp_fac,...
        disp_fac*sim_out_bec_static.final.vel.std(ii)/sqrt(size(sim_out_bec_static.final.vel.vals,1)),'type','b');
end
fprintf('mean vel (x,y,z)= ( %s, %s, %s ) µms⁻¹ \n',str_vals{:})
%set(gcf,'Position',[100 100 1000 1000])

%%

vmag_st=vecnorm(sim_out_bec_static.start.vel.vals,2,2);
vmag_fin_trap_off=vecnorm(sim_out_bec_static.final.vel.vals,2,2);
sout_st=smooth_hist(vmag_st,'scale_y_fun',@(x,y) y.*(x.^4),'sigma',0.001);
sout_fin=smooth_hist(vmag_fin_trap_off,'scale_y_fun',@(x,y) y.*(x.^4),'sigma',0.001);
% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
stfig('simple vel hist')
plot(sout_st.bin.centers,sout_st.count_rate.smooth,'k')
hold on
plot(sout_fin.bin.centers,sout_fin.count_rate.smooth,'b')
hold off
%xlim([opts.vmin,opts.vmax])
%xlim([opts.vmin,opts.vmax])


%% now lets try the simulation with the BEC expanding from trap switch off
% sim_opt=[];
% % sample from TF density squared
% sim_opt.den_power=2;
% % lets use a stationary BEC
% sim_opt.osc.omega =tf_details.inputs.omega;
% %lets define the amplitude in velocity then convert to spatial
% sim_opt.osc.amp=col_vec([0,0,0]); 
% sim_opt.osc.phase=col_vec([0,0,0]);
% sim_opt.nsamp=1e5;
% %sim_opt.grav=const.g0;
% sim_opt.grav=0;
% sim_opt.qd=[];
% % bec is properly finished by 50e-3 a larger inital velociy allows a smaller t_end
% %sim_opt.qd.vmin=50e-3;
% sim_opt.qd.vmin=40e-3;
% sim_opt.qd.vmax=400e-3;
% sim_opt.trap_off=0;
% sim_opt.t_end=0.02;

% use the same simulation options from before but turn the trap off
sim_opt.trap_off=0;

sim_out_trap_off=class_sim_outcoupling(tf_details,sim_opt);
%



%% plots
opts=[];
opts.cords_to_plot=[2,3];
opts.cords_labels={'x','y','z'};
opts.cords_units=repmat({'m $\mathrm{s}^{-1}$'},[1,3]);
opts.mean_marker=true;
opts.disp_scale=[1,1,1]*1e-3;
opts.blur_size=[1,1,1]*0.5e-3;
% opts.std_norm_bounds=true;
% opts.hist_bounds=[5,1,1]*3;
opts.std_norm_bounds=false;
opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
opts.hist_bounds=[[-100,100];[-100,100];[-100,100]]*1e-3;
% optional
opts.filter_in_other_axis=false;
opts.norm_den_to_num_in=true;
opts.scale_den_disp_fac=true;
opts.norm_den_unity=true;
opts.save_hist_as_image=true;
opts.save_hist_name='many_counts.tiff';

%opts.cbar_lims=[0,2e-3];
% opts.bin_factor;

%hist_2d_data(sim_out.start.vel.vals,opts)
hist_2d_data(sim_out_trap_off.final.vel.vals,opts)

hold on
plot(sim_out_trap_off.bec.vel(opts.cords_to_plot(1))/opts.disp_scale(opts.cords_to_plot(1)),...
    sim_out_trap_off.bec.vel(opts.cords_to_plot(2))/opts.disp_scale(opts.cords_to_plot(1)),...
    'wo','MarkerSize',10,'LineWidth',1.0)
hold off

% comapre the mean postion in velocity space to the bec
% this tells us how much the BEC gets pushed arround from the outcoupling process
% 70(7)e-6 ms⁻¹ is what i got in a 1e7 count sim
% the impulse on the bec is therefore 70e-6*nout*mass
% the change in velocity of the bec is therefor 70e-6*nout*mass=m*deltav*nin
% 70e-6 ms⁻¹*nout/nin=deltav
% 70e-6  ms⁻¹*(N_0*eta)/(N_0*(1-eta))=deltav
% for small outcoupling we aproximate this as 
% 70e-6  ms⁻¹*eta=deltav
% for a eta=1e-2
% \delta v=0.7 µs⁻¹


str_vals={};
disp_fac=1e6;
for ii=1:3
    str_vals{ii}=string_value_with_unc(sim_out_trap_off.final.vel.mean(ii)*disp_fac,...
        disp_fac*sim_out_trap_off.final.vel.std(ii)/sqrt(size(sim_out_trap_off.final.vel.vals,1)),'type','b');
end
fprintf('mean vel (x,y,z)= ( %s, %s, %s ) µms⁻¹ \n',str_vals{:})
%set(gcf,'Position',[100 100 1000 1000])

%%


stfig('k log dist')
clf
vmag_st_trap_on=vecnorm(sim_out_bec_static.start.vel.vals,2,2);
vmag_st_trap_off=vecnorm(sim_out_trap_off.start.vel.vals,2,2);
vmag_fin_trap_on=vecnorm(sim_out_bec_static.final.vel.vals,2,2);
vmag_fin_trap_off=vecnorm(sim_out_trap_off.final.vel.vals,2,2);
sigma=1e-3;
sout_st_trap_on=smooth_hist(vmag_st_trap_on,'sigma',sigma); %,'sigma',10
sout_st_trap_off=smooth_hist(vmag_st_trap_off,'sigma',sigma); %,'sigma',10
sout_fin_trap_off=smooth_hist(vmag_fin_trap_off,'sigma',sigma); %v,'sigma',10
sout_fin_trap_on=smooth_hist(vmag_fin_trap_on,'sigma',sigma); %v,'sigma',10

plot(sout_st_trap_on.bin.centers,sout_st_trap_on.count_rate.smooth_prob,'k')
hold on
plot(sout_st_trap_off.bin.centers,sout_st_trap_off.count_rate.smooth_prob,'g')
plot(sout_fin_trap_on.bin.centers,sout_fin_trap_on.count_rate.smooth_prob,'r')
plot(sout_fin_trap_off.bin.centers,sout_fin_trap_off.count_rate.smooth_prob,'b')
hold off
xlabel('$v (\mathrm{m}\cdot\mathrm{s}^{-1})$')
ylabel('density $(\mathrm{m}^{-1} \mathrm{s})$')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
xlim([sim_opt.qd.vmin,sim_opt.qd.vmax])
legend({'start trap on','start trap off','end trap on','end trap off'})


stfig('k neg power mag dist')
clf


% expect this to be -3 for pure QD
scale_power=-(4.0-1);
sigma=10;
% but it seems that -4 is a better fit
scale_power=-4;
sigma=100;

% sout_st=smooth_hist(vmag_st,'scale_y_fun',@(x,y) y.*(x.^4),'sigma',0.001);
% sout_fin=smooth_hist(vmag_fin,'scale_y_fun',@(x,y) y.*(x.^4),'sigma',0.001);
x_scale_fun=@(x) (x.^scale_power);

sout_st_trap_on=smooth_hist(vmag_st_trap_on,'sigma',sigma,'scale_x_fun',x_scale_fun); %,'sigma',10
sout_st_trap_off=smooth_hist(vmag_st_trap_off,'sigma',sigma,'scale_x_fun',x_scale_fun); %,'sigma',10
sout_fin_trap_off=smooth_hist(vmag_fin_trap_off,'sigma',sigma,'scale_x_fun',x_scale_fun); %v,'sigma',10
sout_fin_trap_on=smooth_hist(vmag_fin_trap_on,'sigma',sigma,'scale_x_fun',x_scale_fun); %v,'sigma',10

% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
plot(sout_st_trap_on.bin.centers,sout_st_trap_on.count_rate.smooth_prob,'k')
xl=xlim;
hold on
plot(sout_st_trap_off.bin.centers,sout_st_trap_off.count_rate.smooth_prob,'g')
plot(sout_fin_trap_on.bin.centers,sout_fin_trap_on.count_rate.smooth_prob,'r')
plot(sout_fin_trap_off.bin.centers,sout_fin_trap_off.count_rate.smooth_prob,'b')
hold off
xlim(xl)
xlabel(sprintf('$v^{%.1g}(\\mathrm{m} \\mathrm{s}^{-1})^{%.1g}$',scale_power,scale_power))
ylabel('density')
legend({'start trap on','start trap off','end trap on','end trap off'})

%xlim([opts.vmin,opts.vmax])
%xlim([opts.vmin,opts.vmax])
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')

stfig('k mag dist')
clf
scale_power=4;
y_scale_fun=@(x,y) y.*(x.^scale_power);
sigma=0.001
sout_st_trap_on=smooth_hist(vmag_st_trap_on,'sigma',sigma,'scale_y_fun',y_scale_fun); %,'sigma',10
sout_st_trap_off=smooth_hist(vmag_st_trap_off,'sigma',sigma,'scale_y_fun',y_scale_fun); %,'sigma',10
sout_fin_trap_off=smooth_hist(vmag_fin_trap_off,'sigma',sigma,'scale_y_fun',y_scale_fun); %v,'sigma',10
sout_fin_trap_on=smooth_hist(vmag_fin_trap_on,'sigma',sigma,'scale_y_fun',y_scale_fun); %v,'sigma',10


% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
plot(sout_st_trap_on.bin.centers,sout_st_trap_on.count_rate.smooth_prob,'k')
xl=xlim;
hold on
plot(sout_st_trap_off.bin.centers,sout_st_trap_off.count_rate.smooth_prob,'g')
plot(sout_fin_trap_on.bin.centers,sout_fin_trap_on.count_rate.smooth_prob,'r')
plot(sout_fin_trap_off.bin.centers,sout_fin_trap_off.count_rate.smooth_prob,'b')
hold off
xlabel('$v(\mathrm{m}\cdot\mathrm{s}^{-1})$')
ylabel(sprintf('density $\\cdot k^{%.1g}$',scale_power))
legend({'start trap on','start trap off','end trap on','end trap off'})
%xlim([opts.vmin,opts.vmax])
%xlim([opts.vmin,opts.vmax])

%% lets try filtering out such that we only allow a certian (v_x^2+v_y^2)

mask_rad=0.08;
circle_mask=[0,0,mask_rad,1];
v_filt_fin=sim_out.final.vel.vals;
v_filt_fin(:,:)=v_filt_fin(:,[3,1,2]);% rpermute x,z for the mask function
v_filt_fin=masktxy_2d_circle(v_filt_fin,circle_mask);
v_filt_fin(:,[3,1,2])=v_filt_fin(:,:);% rpermute x,z for the mask function

v_filt_start=sim_out.start.vel.vals;
v_filt_start(:,:)=v_filt_start(:,[3,1,2]);% rpermute x,z for the mask function
v_filt_start=masktxy_2d_circle(v_filt_start,circle_mask);
v_filt_start(:,[3,1,2])=v_filt_start(:,:);% rpermute x,z for the mask function



%%

%% plots
opts=[];
opts.cords_to_plot=[1,2];
opts.cords_labels={'x','y','z'};
opts.cords_units=repmat({'m $\mathrm{s}^{-1}$'},[1,3]);
opts.mean_marker=true;
opts.disp_scale=[1,1,1]*1e-3;
opts.blur_size=[1,1,1]*0.5e-3;
% opts.std_norm_bounds=true;
% opts.hist_bounds=[5,1,1]*3;
opts.std_norm_bounds=false;
opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
opts.hist_bounds=[[-100,100];[-100,100];[-100,100]]*1e-3;
% optional
opts.filter_in_other_axis=false;
opts.norm_den_to_num_in=true;
opts.scale_den_disp_fac=true;
opts.norm_den_unity=true;
opts.save_hist_as_image=true;
opts.save_hist_name='many_counts.tiff';

%opts.cbar_lims=[0,2e-3];
% opts.bin_factor;

%hist_2d_data(sim_out.start.vel.vals,opts)
hist_2d_data(v_filt_start,opts)

hold on
plot(sim_out.bec.vel(opts.cords_to_plot(1))/opts.disp_scale(opts.cords_to_plot(1)),...
    sim_out.bec.vel(opts.cords_to_plot(2))/opts.disp_scale(opts.cords_to_plot(1)),...
    'wo','MarkerSize',10,'LineWidth',1.0)
hold off

% comapre the mean postion in velocity space to the bec
% this tells us how much the BEC gets pushed arround from the outcoupling process
% 70(7)e-6 ms⁻¹ is what i got in a 1e7 count sim
% the impulse on the bec is therefore 70e-6*nout*mass
% the change in velocity of the bec is therefor 70e-6*nout*mass=m*deltav*nin
% 70e-6 ms⁻¹*nout/nin=deltav
% 70e-6  ms⁻¹*(N_0*eta)/(N_0*(1-eta))=deltav
% for small outcoupling we aproximate this as 
% 70e-6  ms⁻¹*eta=deltav
% for a eta=1e-2
% \delta v=0.7 µs⁻¹


str_vals={};
disp_fac=1e6;
for ii=1:3
    str_vals{ii}=string_value_with_unc(sim_out.final.vel.mean(ii)*disp_fac,...
        disp_fac*sim_out.final.vel.std(ii)/sqrt(size(sim_out.final.vel.vals,1)),'type','b');
end
fprintf('mean vel (x,y,z)= ( %s, %s, %s ) µms⁻¹ \n',str_vals{:})
%set(gcf,'Position',[100 100 1000 1000])














%%


stfig('k log dist')
clf
vmag_st=vecnorm(v_filt_start,2,2);
vmag_fin_trap_off=vecnorm(v_filt_fin,2,2);

sout_st=smooth_hist(vmag_st); %,'sigma',10
sout_fin=smooth_hist(vmag_fin_trap_off); %v,'sigma',10

plot(sout_st.bin.centers,sout_st.count_rate.smooth,'k')
hold on
plot(sout_fin.bin.centers,sout_fin.count_rate.smooth,'b')
hold off
xlabel('$v (\mathrm{m}\cdot\mathrm{s}^{-1})$')
ylabel('density $(\mathrm{m}^{-1} \mathrm{s})$')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
xlim([sim_opt.qd.vmin,sim_opt.qd.vmax])


stfig('k neg power mag dist')
clf
vmag_st=vecnorm(sim_out.start.vel.vals,2,2);
vmag_fin_trap_off=vecnorm(sim_out.final.vel.vals,2,2);
% expect this to be -3 for pure QD
scale_power=-(4.0-1);
% but it seems that -4 is a better fit
scale_power=-3;

% sout_st=smooth_hist(vmag_st,'scale_y_fun',@(x,y) y.*(x.^4),'sigma',0.001);
% sout_fin=smooth_hist(vmag_fin,'scale_y_fun',@(x,y) y.*(x.^4),'sigma',0.001);
x_scale_fun=@(x) (x.^scale_power);
sigma=30;
sout_st=smooth_hist(vmag_st,'scale_x_fun',x_scale_fun,'sigma',sigma); %
sout_fin=smooth_hist(vmag_fin_trap_off,'scale_x_fun',x_scale_fun,'sigma',sigma); %v,'sigma',10

% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
plot(sout_st.bin.centers,sout_st.count_rate.smooth,'k')
xl=xlim;
hold on
plot(sout_fin.bin.centers,sout_fin.count_rate.smooth,'b')
hold off
xlim(xl)
xlabel(sprintf('$v^{%.1g}(\\mathrm{m} \\mathrm{s}^{-1})^{%.1g}$',scale_power,scale_power))
ylabel('density')
%xlim([opts.vmin,opts.vmax])
%xlim([opts.vmin,opts.vmax])
% set(gca,'Xscale','log')
% set(gca,'Yscale','log')

stfig('k mag dist')
clf
vmag_st=vecnorm(sim_out.start.vel.vals,2,2);
vmag_fin_trap_off=vecnorm(sim_out.final.vel.vals,2,2);
scale_power=4;
sout_st=smooth_hist(vmag_st,'scale_y_fun',@(x,y) y.*(x.^scale_power),'sigma',0.001);
sout_fin=smooth_hist(vmag_fin_trap_off,'scale_y_fun',@(x,y) y.*(x.^scale_power),'sigma',0.001);


% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
plot(sout_st.bin.centers,sout_st.count_rate.smooth,'k')
hold on
plot(sout_fin.bin.centers,sout_fin.count_rate.smooth,'b')
hold off
xlabel('$v(\mathrm{m}\cdot\mathrm{s}^{-1})$')
ylabel(sprintf('density $\\cdot k^{%.1g}$',scale_power))

%xlim([opts.vmin,opts.vmax])
%xlim([opts.vmin,opts.vmax])

