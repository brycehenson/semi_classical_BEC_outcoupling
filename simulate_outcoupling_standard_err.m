%simulate standard error in fit

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path
hebec_constants
rng('shuffle')


%% simulate the standar error in the fit frequency by repeating many times
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
osc_amp_vel=col_vec([0,1,0]*10e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.lambda=0.30;%5/500;
sim_opt.osc.phase=col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.grav=const.g0;
sim_opt.nsamp=3.2e3*0.1;%number of atoms to simulate


%times=linspace(0,50,156).*(2*pi/sim_opt.osc.omega(2));
times=(0:(156-1))*8e-3; %156
jjmax=numel(times);


iimax=10;
fit_repeats=[];
fit_repeats.fit.freq.val=[];
fit_repeats.fit.freq.se=[];
fit_repeats.fitparam={};


%sim_opt.verbose=3;
%sim_outs=cell(1,jjmax);

for ii=1:iimax

    vel_compare=[];
    vel_compare.bec.mean=nan(jjmax,3);
    vel_compare.pal.mean=nan(jjmax,3);
    vel_compare.pal.std=nan(jjmax,3);
    vel_compare.pal.ste=nan(jjmax,3);
    vel_compare.time=times;

    for jj=1:jjmax
        % do the simulation of a out-coupling time
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
    end
    
    
    idx=2;
    data_in=[];
    data_in.x=vel_compare.pal.mean(:,idx);
    data_in.t=vel_compare.time;
    data_in.unc_x=vel_compare.pal.ste(:,idx)
    opts_in=[];
    fit_pal=fit_sine_to_data(data_in,opts_in)

    fit_repeats.fit.freq.val(ii)=fit_pal.fitparam.Estimate(3);
    fit_repeats.fit.freq.se=fit_pal.fitparam.SE(3);
    fit_repeats.fitparam{ii}=fit_pal.fitparam;

  
end

%%
std(rmoutliers(fit_repeats.fit.freq.val))
std_err_of_sample_std(numel(rmoutliers(fit_repeats.fit.freq.val)),std(rmoutliers(fit_repeats.fit.freq.val)))


%%
stfig('pal as momentum sampler')
y_labels_val={'V x (mm/s)','V y  (mm/s)','V z  (mm/s)'}
y_labels_diff={'$\Delta$ V x (mm/s)','$\Delta$ V y  (mm/s)','$\Delta$ V z  (mm/s)'}
font_size=15
y_factor=1e3
for idx=1:3
    subplot(3,2,2*(idx-1)+1)
    plot(vel_compare.time,vel_compare.bec.mean(:,idx)*y_factor,'k')
    ylabel(y_labels_val{idx},"FontSize",font_size)
    xlabel('Time (s)',"FontSize",font_size)
    hold on
    plot(vel_compare.time,vel_compare.pal.mean(:,idx)*y_factor)
    hold off
    legend('a','b')
    subplot(3,2,2*(idx-1)+2)
    plot(vel_compare.time,(vel_compare.bec.mean(:,idx)-vel_compare.pal.mean(:,idx)) *y_factor )
    ylabel(y_labels_diff{idx},"FontSize",font_size)
    xlabel('Time (s)',"FontSize",font_size)
end
ax=gca
