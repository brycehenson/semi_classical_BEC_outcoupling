%simulate_outcoupling

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path

hebec_constants
rng('shuffle')

%%
num_local=change_num_workers_local()

parpool('local',num_local)

%%




%% lets study what happens when there is an oscillation of the condensate
% we will usea a 10mm/s oscillation because thats about what the tune out used

%define exp parameters
tf_param=[];
tf_param.trap_freq=2*pi*[55,426,428];
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
% lets define a time varying atom number
tf_num_tume_fun=@(t) 6e5-(5e5/1.2)*t;

do_plot=false;
sim_opt=[];
sim_opt.osc.omega =tf_details.inputs.omega;
%lets define the amplitude in velocity then convert to spatial
osc_amp_vel=col_vec([0,1,0]*12e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.lambda=1/0.7;%5/500;
sim_opt.osc.phase=col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.grav=const.g0;
sim_opt.nsamp=3.2e3*0.1;%number of atoms to simulate


%times=linspace(0,50,156).*(2*pi/sim_opt.osc.omega(2));
times=(0:(156-1))*8e-3; %156
jjmax=numel(times);

vel_compare=[];
vel_compare.bec.mean=nan(jjmax,3);
vel_compare.pal.mean=nan(jjmax,3);
vel_compare.pal.std=nan(jjmax,3);
vel_compare.pal.ste=nan(jjmax,3);
vel_compare.time=times;


%sim_opt.verbose=3;
%sim_outs=cell(1,jjmax);

for jj=1:jjmax
fprintf('\n time %u of %u \n',jj,jjmax) 

tf_param.num=tf_num_tume_fun(times(jj));
tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);
sim_opt.osc.time_offset=times(jj);
sim_out=class_sim_outcoupling(tf_details,sim_opt);
%sim_outs{jj}=sim_out;
 
vel_compare.bec.mean(jj,:)=row_vec(sim_out.bec.vel);
vel_compare.pal.mean(jj,:)=sim_out.final.vel.mean;
vel_compare.pal.std(jj,:)=sim_out.final.vel.std;
vel_compare.pal.ste(jj,:)=sim_out.final.vel.std/...
                            sqrt(size(sim_out.final.vel.vals,1));

% lets fit the motion that was simulated



if do_plot
    % plots
    hist_opt=[];
    hist_opt.cords_to_plot=[2,3];
    hist_opt.cords_labels={'x','y','z'};
    hist_opt.cords_units=repmat({'m $\mathrm{s}^{-1}$'},[1,3]);
    hist_opt.mean_marker=true;
    hist_opt.disp_scale=[1,1,1]*1e-3;
    hist_opt.blur_size=[1,1,1]*0.2e-3;
    % opts.std_norm_bounds=true;
    % opts.hist_bounds=[5,1,1]*3;
    hist_opt.std_norm_bounds=false;
    hist_opt.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
    hist_opt.hist_bounds=[[-50,50];[-55,55];[-50,100]]*1e-3;
    % optional
    hist_opt.filter_in_other_axis=false;
    hist_opt.norm_den_to_num_in=true;
    hist_opt.scale_den_disp_fac=true;
    %opts.cbar_lims=[0,2e-3];
    % opts.bin_factor;

    hist_2d_data(sim_out.final.vel.vals,hist_opt)

    hold on
    plot(sim_out.bec.vel(hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
        sim_out.bec.vel(hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
        'wo','MarkerSize',10,'LineWidth',1.0)
    hold off
    pause(1e-6)
end
% set(gcf,'Position',[100 100 1000 1000])
% saveas(gcf,['./figs/phase_seq/',sprintf('10_mm_s_phase_pi_%.3f.png',phases(jj))])

end

%save('done_simulating')

%%
stfig('pal as momentum sampler')
y_labels_val={'V x (mm/s)','V y  (mm/s)','V z  (mm/s)'}
y_labels_diff={'$\Delta$ V x (mm/s)','$\Delta$ V y  (mm/s)','$\Delta$ V z  (mm/s)'}
font_size=15
y_factor=1e3
for ii=1:3
    subplot(3,2,2*(ii-1)+1)
    plot(vel_compare.time,vel_compare.bec.mean(:,ii)*y_factor,'k')
    ylabel(y_labels_val{ii},"FontSize",font_size)
    xlabel('Time (s)',"FontSize",font_size)
    hold on
    plot(vel_compare.time,vel_compare.pal.mean(:,ii)*y_factor)
    hold off
    legend('a','b')
    subplot(3,2,2*(ii-1)+2)
    plot(vel_compare.time,(vel_compare.bec.mean(:,ii)-vel_compare.pal.mean(:,ii)) *y_factor )
    ylabel(y_labels_diff{ii},"FontSize",font_size)
    xlabel('Time (s)',"FontSize",font_size)
end
ax=gca

%%
tf_param=[];
tf_param.trap_freq=2*pi*[55,426,428];
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
tf_param.num=4e5;

tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);

unc_est_st=[];
unc_est_st.pal_num_per_pulse=3.2e3*0.1;
unc_est_st.num_pulses=156;
unc_est_st.osc_amp=12e-3;
unc_est_st.time=156*8e-3;
unc_est_st.damp_rate=1/0.7;
out=unc_est_in_pal_meas(tf_details,unc_est_st)


%%
ii=2
data_in=[];
data_in.x=vel_compare.pal.mean(:,ii);
data_in.t=vel_compare.time;
data_in.unc_x=vel_compare.pal.ste(:,ii)
opts_in=[];
fit_pal=fit_sine_to_data(data_in,opts_in)

%%

data_in=[];
data_in.x=vel_compare.bec.mean(:,ii);
data_in.t=vel_compare.time;
opts_in=[];
fit_bec=fit_sine_to_data(data_in,opts_in)


%% look for harnomics
% we will repeat the data a number of times
repeats=100;
dim_idx=2;
osc_period=1/(sim_opt.osc.omega(dim_idx)/(2*pi));
time_vel=zeros(0,2);
frac_cyc=vel_compare.phase/(2*pi);
frac_cyc=col_vec(frac_cyc);
pal_vel=col_vec(vel_compare.pal(:,dim_idx));
for ii=0:(repeats-1)
    tmp_time_vel=cat(2,osc_period*(frac_cyc+ii),pal_vel);
    time_vel=cat(1,time_vel,tmp_time_vel);
end

stfig('repeated motion')
plot(time_vel(:,1),time_vel(:,2))
% if we have done the shift right there should be no difference if we shift the data an integer num of periods
hold on
plot(time_vel(:,1)+3*osc_period,time_vel(:,2))
hold off

%%
fft_out=fft_tx(time_vel(:,1),time_vel(:,2)-mean(time_vel(:,2)),'padding',3,'window','chebyshev','win_param',{100});
%out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
fft_amp=abs(fft_out(2,:));
max_fft_amp=max(fft_amp);
plot(fft_out(1,:)*1e-3,1e3*fft_amp)
set(gca,'Yscale','log')
xlabel('freq (kHz)')
ylabel('amplitude (mm)')
xlim([0,5.2])
ylim([1e-4,15])

%%
hist_opt=[];
hist_opt.cords_to_plot=[2,3];
hist_opt.cords_labels={'x','y','z'};
hist_opt.blur_size=[1,1,1]*1e3;
%outcoupled.vel
hist_spatial_data(outcoupled.pos,hist_opt)



%%


%%
[jac,jac_err] = jacobianest(@(x) state_update_fun_grav(x,tf_details,xshape,const.g0) ,xstart)

d_dxdt_dx=state_update_jacobian(xstart,xshape,tf_details)

%%
state_update_fun_grav(xstart,tf_details,xshape,const.g0)

%%


