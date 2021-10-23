%simulate_outcoupling

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path

hebec_constants
rng('shuffle')

%%
num_local=change_num_workers_local(15)

parpool('local',num_local)
%%

%define exp parameters
param=[];
param.trap_freq=2*pi*[50,500,500];
param.num=1e6;
param.mass=const.mhe;
param.a_scat_len=const.ahe_scat;

tf_details=bec_properties(param.trap_freq/(2*pi),param.num,param.mass,param.a_scat_len);


%% lets study what happens when there is an oscillation of the condensate
% we will usea a 10mm/s oscillation because thats about what the tune out used
do_plot=true;
sim_opt=[];
sim_opt.osc.omega =tf_details.inputs.omega;
%lets define the amplitude in velocity then convert to spatial
osc_amp_vel=col_vec([0,1,1]*12e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.decay_tau=inf;%5/500;
sim_opt.osc.phase=col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.nsamp=3e6;%number of atoms to simulate

%phases=linspace(0,2,10)*pi; %2^7 2^12
%phases=phases(1:end-1); % remove the 2 pi pt because its equal to the start
%jjmax=numel(phases);

times=linspace(0,1,200).*(2*pi/sim_opt.osc.omega(2));
times(end)=[]; % for a nice loop remove the last point
jjmax=numel(times);


vel_compare=[];
vel_compare.bec.mean=nan(jjmax,3);
vel_compare.pal.mean=nan(jjmax,3);
vel_compare.pal.std=nan(jjmax,3);
vel_compare.pal.ste=nan(jjmax,3);
vel_compare.time=times;


save_path='E:/sratch/pretty_3/';

%sim_opt.verbose=3;
%sim_outs=cell(1,jjmax);
jjmin=1;
for jj=jjmin:jjmax
    fprintf('\n time %u of %u \n',jj,jjmax) 
    %sim_opt.osc.phase=col_vec([0,0,(1/2)])*pi+ phases(jj);
    sim_opt.osc.time_offset=times(jj);
    
    sim_out={}; % free up memory
    sim_out=class_sim_outcoupling(tf_details,sim_opt);
    
    %sim_outs{jj}=sim_out;
     
    vel_compare.bec.mean(jj,:)=row_vec(sim_out.bec.vel);
vel_compare.pal.mean(jj,:)=sim_out.final.vel.mean;
vel_compare.pal.std(jj,:)=sim_out.final.vel.std;
vel_compare.pal.ste(jj,:)=sim_out.final.vel.std/...
                            sqrt(size(sim_out.final.vel.vals,1));
    
    if do_plot
        % plots
        hist_opt=[];
        hist_opt.cords_to_plot=[2,3];
        hist_opt.cords_labels={'x','y','z'};
        hist_opt.cords_units=repmat({'m $\mathrm{s}^{-1}$'},[1,3]);
        hist_opt.mean_marker=true;
        hist_opt.disp_scale=[1,1,1]*1e-3;
        hist_opt.blur_size=[1,1,1]*0.3e-3; %0.3e-3
        hist_opt.kernel_factor=10;
        % opts.std_norm_bounds=true;
        % opts.hist_bounds=[5,1,1]*3;
        hist_opt.std_norm_bounds=false;
        %opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
        hist_opt.cmap=viridis(512);
        hist_opt.density_power=0.4;
        hist_opt.hist_bounds=[[-50,50];[-60,60];[-50,100]]*1e-3;
        % optional
        hist_opt.filter_in_other_axis=false;
        hist_opt.norm_den_to_num_in=true;
        hist_opt.scale_den_disp_fac=true;
        %opts.cbar_lims=[0,2e-3];
        % opts.bin_factor;
        hist_opt.save_hist_as_image=true;
        hist_opt.save_hist_name=strcat(save_path,sprintf('raw_%0.4u.png',jj));
    
    
        hist_fig=hist_2d_data(sim_out.final.vel.vals,hist_opt)
        
        hold on
        plot(sim_out.bec.vel(hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
            sim_out.bec.vel(hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
            'wo','MarkerSize',10,'LineWidth',1.0)
        hold off
        pause(1e-6)
    end
    
    set(gcf,'Position',[100 100 1000 1000])
    if ~strcmp(save_path,'')
        saveas(hist_fig,[save_path,sprintf('anotated_step_%0.4u.png',jj)])
    end

end

%save('done_simulating')

%% make into video
% cut some music down
% ffmpeg -ss 3:40 -t 20 -i music.mp3 short_music.mp3 

% new workflow
% ffmpeg -f image2 -framerate 12 -i raw_%04d.png -filter:v "minterpolate='mi_mode=blend:fps=24',spp" -vcodec libx264 -b 16000k  single_seq.mp4
% loop ffmpeg -stream_loop 1 -i single_seq.mp4 -c copy loop_seq.mp4 


% add it to video
% ffmpeg -f image2 -framerate 12 -i raw_%04d.png -i short_music.mp3 -shortest -filter:v "minterpolate='mi_mode=blend:fps=24',spp" -vcodec libx264 -b 16000k  video.mp4
%%
stfig('pal as momentum sampler')
for ii=1:3
subplot(3,2,2*(ii-1)+1)
plot(vel_compare.phase,vel_compare.bec(:,ii))
hold on
plot(vel_compare.phase,vel_compare.pal(:,ii))
hold off
subplot(3,2,2*(ii-1)+2)
plot(vel_compare.phase,vel_compare.bec(:,ii)-vel_compare.pal(:,ii))
end


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


