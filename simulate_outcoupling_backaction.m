% find the impulse which is imparted on the BEC from outcoupling

addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path
hebec_constants
rng('shuffle')

%%
num_local=change_num_workers_local(20)


%% Measure the impulse on the BEC from outcoupling



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
%osc_amp_vel=col_vec([0,0,0.00]*1e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.lambda=0.30;%5/500;
sim_opt.osc.phase=col_vec([1/4,0,(0)])*pi; %col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.grav=const.g0;
sim_opt.nsamp=8e5;%number of atoms to simulate
sim_opt.osc.time_offset=0;
sim_opt.samp_dynamics_times=linspace(0,8e-3,200);
sim_opt.t_end=12e-3;
sim_opt.verbose=2;


% do outoupling simulation
sim_out=class_sim_outcoupling(tf_details,sim_opt);

 
%vel_compare.bec.mean(jj,:)=row_vec(sim_out.bec.vel);
%vel_compare.pal.mean(jj,:)=sim_out.final.vel.mean;
%vel_compare.pal.std(jj,:)=sim_out.final.vel.std;
%vel_compare.pal.ste(jj,:)=sim_out.final.vel.std/...
%                            sqrt(size(sim_out.final.vel.vals,1));

%%

deltav= sim_out.start.vel.mean-sim_out.final.vel.mean;
impulse=deltav*const.mhe;

signifigance=deltav./sim_out.final.vel.ste;


%% velocity



font_name='cmr10';
linewidth=1.5;
font_size=12;

colors_main=[[210,56,46];[50,163,73];[73,133,173]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);

y_factor=1e3;
x_factor=1e3;
%x_lim=[min(vel_compare.time*x_factor),max(vel_compare.time*x_factor)];
x_lim=[-0.1,8];

stfig('mean vel')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
subplot(1,2,1)
shaded_ci_lines=true;
hold on
line_style={':','-','--'}
time_scale=1e3;
yscale=1e3;
for idx=1:3
    plot(sim_out.dyn_states.times*time_scale,(sim_out.dyn_states.vel.mean(:,idx)-sim_out.bec.vel(idx))*yscale,...
        line_style{idx},'LineWidth',linewidth,'Color',colors_main(idx,:))
    if shaded_ci_lines
        ci_up=(sim_out.dyn_states.vel.mean(:,idx)-sim_out.bec.vel(idx)+sim_out.dyn_states.vel.ste(:,idx))*yscale;
        ci_down=(sim_out.dyn_states.vel.mean(:,idx)-sim_out.bec.vel(idx)-sim_out.dyn_states.vel.ste(:,idx))*yscale;
        patch([sim_out.dyn_states.times*time_scale', fliplr(sim_out.dyn_states.times*time_scale')], ...
            [ci_up', fliplr(ci_down')], colors_shaded(idx,:),'EdgeColor','none',...
            'FaceAlpha',0.3)  %[1,1,1]*0.80
    else
        plot(sim_out.dyn_states.times*time_scale,(sim_out.dyn_states.vel.mean(:,idx)-sim_out.dyn_states.vel.ste(:,idx))*yscale,...
            line_style{idx},'LineWidth',1.0,'Color',colors_shaded(idx,:))
        plot(sim_out.dyn_states.times*time_scale,(sim_out.dyn_states.vel.mean(:,idx)+sim_out.dyn_states.vel.ste(:,idx))*yscale,...
            line_style{idx},'LineWidth',1.0,'Color',colors_shaded(idx,:))
    end
end
if shaded_ci_lines
    legend('x','','y','','z','Location','best','FontSize',font_size*1.25)
else
    legend('x','','','y','','','z','Location','best')
end
legend('boxoff')
legend('FontSize',font_size*1.25)
hold off
box on
xlabel('time (ms)')
ylabel('mean velocity (mm$\cdot\mathrm{s}^{-1}$)')
set(gca,'xlim',x_lim)
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])

sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.01, sp_pos(2)+sp_pos(4)*0.99, 0, 0], 'string', '(A)','FontSize',round(font_size*1.2),'FontName',font_name);


% acceleration

% calculate the derivative

sim_out.dyn_states.acc.mean=diff(sim_out.dyn_states.vel.mean,1,1)./repmat(diff(sim_out.dyn_states.times,1,2)',[1,3]);
sim_out.dyn_states.acc.ste=sqrt(sim_out.dyn_states.times(2:end).^2+sim_out.dyn_states.times(1:end-1).^2)./diff(sim_out.dyn_states.times,1,2);
sim_out.dyn_states.acc.time=(1/2)*(sim_out.dyn_states.times(2:end)+sim_out.dyn_states.times(1:end-1));


subplot(1,2,2)
hold on
time_scale=1e3;
yscale=1;
for idx=1:3
    plot(sim_out.dyn_states.acc.time*time_scale,sim_out.dyn_states.acc.mean(:,idx)*yscale,...
        line_style{idx},'LineWidth',linewidth,'Color',colors_main(idx,:))
     ci_up=(sim_out.dyn_states.acc.mean(:,idx)+sim_out.dyn_states.acc.ste(:,idx))*yscale;
     ci_down=(sim_out.dyn_states.acc.mean(:,idx)-sim_out.dyn_states.acc.ste(:,idx))*yscale;
%     if shaded_ci_lines
%         patch([sim_out.dyn_states.acc.time*time_scale', fliplr(sim_out.dyn_states.acc.time*time_scale')], ...
%             [ci_up', fliplr(ci_down')], color_shaded(idx,:),'EdgeColor','none',...
%             'FaceAlpha',0.3)  %[1,1,1]*0.80
%     else
%         plot(sim_out.dyn_states.acc.time*time_scale,ci_up,...
%             line_style{idx},'LineWidth',1.0,'Color',colors_detail(idx,:))
%         plot(sim_out.dyn_states.acc.time*time_scale,ci_down,...
%             line_style{idx},'LineWidth',1.0,'Color',colors_detail(idx,:))
%     end
end
% if shaded_ci_lines
%     legend('x','','y','','z','Location','best','FontSize',15)
% else
%     legend('x','','','y','','','z','Location','southeast')
% end
legend('x','y','z','Location','southeast','FontSize',font_size*1.25)
legend('boxoff')
hold off
box on
xlabel('time (ms)')
ylabel('mean acceleration (m$\cdot\mathrm{s}^{-2}$)')
set(gca,'xlim',x_lim)


set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=1000;
fig_aspect_ratio=0.35; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.01, sp_pos(2)+sp_pos(4)*0.99, 0, 0], 'string', '(B)','FontSize',round(font_size*1.2),'FontName',font_name);


%grid on
fig_name='pal_mean_v_acc_dyn_osc';
%fig_name='pal_mean_v_acc_dyn_static';
%fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
%exportgraphics(gca,fullfile(fig_dir,strcat(fig_name,'.pdf')))




 


%% plot the 2d vel profile of the atom laser in velocity for a number of query points
time_querys=[2.50,4.25,6,12]*1e-3;
time_idxs=nan*time_querys;
for ii=1:numel(time_querys)
    time_idxs(ii)=nth_fun_output(2,@min,abs(sim_out.dyn_states.times-time_querys(ii)));
end 

fh=stfig('hists vel');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

iimax=numel(time_idxs);
for ii=1:iimax
    fh_s=subplot(1,iimax,ii);
    
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
    hist_opt.hist_bounds=[[-40,40];[-40,40];[-30,80]]*1.1*1e-3;
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
    
    
    hist_fig=hist_2d_data( sim_out.dyn_states.vel.vals{time_idxs(ii)},hist_opt);
    
    
    hold on
    plot(sim_out.bec.vel(hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
        sim_out.bec.vel(hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
        '^','MarkerSize',10,'LineWidth',1.4,'Color',[0.1,0.1,1])
    plot(sim_out.dyn_states.bec.vel(time_idxs(ii),hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
        sim_out.dyn_states.bec.vel(time_idxs(ii),hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
        'wo','MarkerSize',5,'LineWidth',1.3)
    hold off

    if ii~=iimax
        ch=colorbar;
        delete(ch)
    end

     string_time=sprintf('t=%.1f ms',sim_out.dyn_states.times(time_idxs(ii))*1e3);
    sp_pos=get(gca,'Position');
    an=annotation('textbox', [sp_pos(1)+0.04, sp_pos(2)+sp_pos(4)*0.95, 0.08, 0.001], 'string', string_time,...
        'FontSize',round(font_size*1.2),'FontName',font_name,'Color',[1,1,1]);
    an.FitBoxToText='on';
    an.LineStyle='none';
end

set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=2000;
fig_aspect_ratio=0.16; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.04,0])
%grid on
fig_name='pal_hist_dyn_vel';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))


%%

%% plot the 2d spatial profile of the atom laser in velocity for a number of query points



fh=stfig('hists pos');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

iimax=numel(time_idxs);
for ii=1:iimax
    fh_s=subplot(1,iimax,ii);
    time_this=sim_out.dyn_states.times(time_idxs(ii));

    % plots
    hist_opt=[];
    hist_opt.cords_to_plot=[2,3];
    hist_opt.cords_labels={'x','y','z'};
    hist_opt.cords_units=repmat({'m'},[1,3]);
    hist_opt.mean_marker=true;
    hist_opt.disp_scale=[1,1,1]*1e-3;
    
    hist_opt.kernel_factor=10;
    % opts.std_norm_bounds=true;
    % opts.hist_bounds=[5,1,1]*3;
    hist_opt.std_norm_bounds=false;
    %opts.cmap=nonlinear_colormap(viridis,'power',[0.45,1,0]);
    hist_opt.cmap=viridis(512);
    hist_opt.density_power=0.4;
    %tf_details.pal_vmax
    hist_opt.hist_bounds=[[-40,40];[-40,40];[-30,80]]*1e-3*time_this;
    
    hist_opt.blur_size=[1,1,1]*2.5e-4*time_this; %0.3e-3
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
    
    
    hist_fig=hist_2d_data( sim_out.dyn_states.pos.vals{time_idxs(ii)},hist_opt);
    
    
    hold on
%     plot(sim_out.bec.vel(hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
%         sim_out.bec.vel(hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
%         'wo','MarkerSize',10,'LineWidth',1.0)
    plot(sim_out.dyn_states.bec.pos(time_idxs(ii),hist_opt.cords_to_plot(1))/hist_opt.disp_scale(hist_opt.cords_to_plot(1)),...
        sim_out.dyn_states.bec.pos(time_idxs(ii),hist_opt.cords_to_plot(2))/hist_opt.disp_scale(hist_opt.cords_to_plot(2)),...
        'wo','MarkerSize',5,'LineWidth',1.3)
    hold off

    if ii~=iimax
        ch=colorbar;
        delete(ch)
    end

    string_time=sprintf('t=%.1f ms',sim_out.dyn_states.times(time_idxs(ii))*1e3);
    sp_pos=get(gca,'Position');
    an=annotation('textbox', [sp_pos(1)+0.04, sp_pos(2)+sp_pos(4)*0.95, 0.08, 0.001], 'string', string_time,...
        'FontSize',round(font_size*1.2),'FontName',font_name,'Color',[1,1,1]);
    an.FitBoxToText='on';
    an.LineStyle='none';
end

set(gca, {'XColor', 'YColor'}, {[1,1,1]*0.2, [1,1,1]*0.2});
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
%fig_width_px=2000;
%fig_aspect_ratio=0.16; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.04,0])
%grid on
fig_name='pal_hist_dyn_pos';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))

