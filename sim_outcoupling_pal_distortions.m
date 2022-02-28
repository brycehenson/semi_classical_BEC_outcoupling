% sim and compare pal distortions



%define exp parameters
% we will usea a 10mm/s oscillation because thats about what the tune out used
tf_param=[];
%tf_param.trap_freq=2*pi*[51,412,415]; 
tf_param.trap_freq=2*pi*[51,446.28,450.7]; 
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
tf_param.num=1.7e5;
tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);

sim_opt=[];
sim_opt.osc.omega =tf_details.inputs.omega;
%lets define the amplitude in velocity then convert to spatial
%osc_amp_vel=col_vec([0,10,0.00]*1e-3); 
osc_amp_vel=col_vec([0,11.1,1.3]*1e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.lambda=0.30;%5/500;
sim_opt.osc.phase=col_vec([0,0.5786,0.658])*pi; %col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.grav=const.g0;
sim_opt.nsamp=2e5;%number of atoms to simulate
sim_opt.osc.time_offset=nan;
sim_opt.t_end=12e-3;
sim_opt.verbose=2;


out_times=(0:3)*7.5e-3;
iimax=numel(out_times);
sim_outs=cell(iimax,1);
for ii=1:iimax
    sim_opt.osc.time_offset=out_times(ii);
    % do outoupling simulation
    sim_outs{ii}=class_sim_outcoupling(tf_details,sim_opt);
end 



%% plot the 2d vel profile of the atom laser at these times

fh=stfig('hists vel');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

iimax=numel(sim_outs);
for ii=1:iimax
    fh_s=subplot(1,iimax,ii);
    sim_out_tmp=sim_outs{ii};
    
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
    hist_opt.hist_bounds=[[-40,40];[-40,40];[-25,47]]*1.1*1e-3;
    hist_opt.blur_size=[1,1,1]*5e-4; %0.3e-3
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
    
    
    hist_fig=hist_2d_data( sim_out_tmp.final.vel.vals,hist_opt);
    
    
    hold on
%     plot(sim_out_tmp.bec.vel(hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
%         sim_out_tmp.bec.vel(hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
%         '^','MarkerSize',10,'LineWidth',1.4,'Color',[0.1,0.1,1])
    plot(sim_out_tmp.bec.vel(hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
        sim_out_tmp.bec.vel(hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
        'wo','MarkerSize',5,'LineWidth',1.3)
    hold off

    if ii~=iimax
        ch=colorbar;
        delete(ch)
    end

    string_time=sprintf('t=%.1f ms',out_times(ii)*1e3);
    sp_pos=get(gca,'Position');
    an=annotation('textbox', [sp_pos(1)+sp_pos(3)*0.2, sp_pos(2)+sp_pos(4)*0.95, 0.08, 0.001], 'string', string_time,...
        'FontSize',round(font_size*1.2),'FontName',font_name,'Color',[1,1,1]);
    an.FitBoxToText='on';
    an.LineStyle='none';
    set(gca,'linewidth', 1.1)
    set(gca,'TickLength',[0.04,0])
    set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});

end

set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=2000;
fig_aspect_ratio=0.15; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])
%grid on
fig_name='pal_hist_vel_out_times_sim';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))