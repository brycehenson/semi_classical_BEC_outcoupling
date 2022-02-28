% simulate outcoupling observed phase shift 

% the outcoupling has a phase shift because of the secondary collision with
% the oscillating bec. make a plot of what this phase shift is as a
% function of atom number.


addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be './lib/Core_BEC_Analysis/lib/'
set_up_project_path
hebec_constants
rng('shuffle')


%% simulate the error in the systematic velocity error as a function of time
% to find the apparent phase shift

%define exp parameters
tf_param=[];
%tf_param.trap_freq=2*pi*[55,426,428];
tf_param.trap_freq=2*pi*[51,412,415];
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
% lets define a time varying atom number
%tf_num_tume_fun=@(t) 6e5-(5e5/1.2)*t;
tf_num_tume_fun=@(t) 6e5;

do_plot=false;
sim_opt=[];
sim_opt.osc.omega =tf_details.inputs.omega;
%lets define the amplitude in velocity then convert to spatial
osc_amp_vel=col_vec([0,10,0]*1e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.lambda=0;% 0.30;%5/500;
sim_opt.osc.phase=col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.grav=const.g0;
sim_opt.nsamp=1.5e4;%number of atoms to simulate


times=linspace(0,3*(2*pi/sim_opt.osc.omega(2)),100);
jjmax=numel(times);

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


fit_idx=2;
data_in=[];
data_in.x=vel_compare.pal.mean(:,fit_idx);
data_in.t=vel_compare.time;
data_in.unc_x=vel_compare.pal.ste(:,fit_idx);
opts_in=[];
fit_pal=fit_sine_to_data(data_in,opts_in);
data_in=[];
data_in.x=vel_compare.bec.mean(:,fit_idx);
data_in.t=vel_compare.time;
fit_bec=fit_sine_to_data(data_in,opts_in);

res_st={};
res_st.phase_shift={};
res_st.phase_shift.val=fit_pal.fitparam{'phase','Estimate'}-fit_bec.fitparam{'phase','Estimate'};
res_st.phase_shift.se=sqrt( fit_pal.fitparam{'phase','SE'}^2+fit_bec.fitparam{'phase','SE'}^2);

rad2deg(res_st.phase_shift.val*2*pi)
fit_pal.fitparam{'amp','Estimate'}
fit_pal.fitparam{'amp','SE'}

%% plot this sim

font_name='cmr10';
linewidth=1.2;
font_size=12;

colors_main=[[210,56,46];[50,163,73];[73,133,173]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);

y_factor=1e3;
x_factor=1e3;
x_lim=[min(vel_compare.time*x_factor),max(vel_compare.time*x_factor)];


stfig('velocity error')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

subplot(2,1,1)
plot_handles={};
hold on
box on
for idx=1:3
    plot_handles{idx}=plot(vel_compare.time*x_factor,vel_compare.pal.mean(:,idx)*y_factor, ...
        'color',colors_main(idx,:),'LineWidth',linewidth);
end
for idx=1:3
    plot_handles{idx+3}=plot(vel_compare.time*x_factor,vel_compare.bec.mean(:,idx)*y_factor,'--',...
        'color',colors_main(idx,:),'LineWidth',linewidth);
end
hold off
xlim(x_lim)
ln=legend([plot_handles{:}] ...
        ,'$x$ PAL','$y$ PAL','$z$ PAL','$x$ BEC','$y$ BEC','$z$ BEC')
legend('FontSize',font_size*0.7)
legend('Location','best')
pause(0.1)
ln.Position=ln.Position
legend('Location','best')
ylabel('velocity, $v$ (mm/s)',"FontSize",font_size)
xlabel('Time (ms)',"FontSize",font_size)
% text(2.2,12,'$y$','FontSize',25)
% text(2.2,5.5,'$z$','FontSize',25)
% text(2.2,-2,'$x$','FontSize',25)
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.1, sp_pos(2)+sp_pos(4)*0.95, 0, 0], 'string', '(A)','FontSize',round(font_size*1.2),'FontName',font_name);


subplot(2,1,2)
box on
hold on
plot_handles={};
line_style={':','-','--'};
% plot the error bars
for idx=1:3
    %plot(vel_compare.time,(vel_compare.pal.mean(:,idx)-vel_compare.bec.mean(:,idx))*y_factor,'k')
    x_val=vel_compare.time';
    y_val=(vel_compare.pal.mean(:,idx)-vel_compare.bec.mean(:,idx));
    err_pos=y_val+vel_compare.pal.ste(:,idx);
    err_neg=y_val-vel_compare.pal.ste(:,idx);
    plot_handles{idx+3}=...
        fill([x_val; flipud(x_val)].*x_factor,[err_neg; flipud(err_pos)].*y_factor,...
            colors_shaded(idx,:),'EdgeAlpha',0);
    set(plot_handles{idx+3},'facealpha',.3)
        
end
for idx=1:3
    %plot(vel_compare.time,(vel_compare.pal.mean(:,idx)-vel_compare.bec.mean(:,idx))*y_factor,'k')
    x_val=vel_compare.time';
    y_val=(vel_compare.pal.mean(:,idx)-vel_compare.bec.mean(:,idx));
    plot_handles{idx}=...
        plot(x_val*x_factor,y_val*y_factor,'color',colors_main(idx,:), ...
        'LineWidth',linewidth,'LineStyle',line_style{idx});
        
end
hold off
xlim(x_lim)
ln=legend([plot_handles{1:3}],'$x$','$y$','$z$')
legend('Location','best')
legend('FontSize',font_size*0.8)
pause(0.1)
ln.Position=ln.Position
ylabel('velocity difference,$\Delta v$ (mm/s)',"FontSize",font_size)
xlabel('Time (ms)',"FontSize",font_size)
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.1, sp_pos(2)+sp_pos(4)*0.95, 0, 0], 'string', '(B)','FontSize',round(font_size*1.2),'FontName',font_name);


set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=700;
fig_aspect_ratio=0.8; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])
%grid on
fig_name='pal_apparent_velocity';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
%exportgraphics(gca,fullfile(fig_dir,strcat(fig_name,'.pdf')))





%% loop over many atom numbers


%define exp parameters
tf_param=[];
%tf_param.trap_freq=2*pi*[55,426,428];
tf_param.trap_freq=2*pi*[51,412,415]';
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;

do_plot=false;
sim_opt=[];
sim_opt.osc.omega =tf_param.trap_freq;
%lets define the amplitude in velocity then convert to spatial
osc_amp_vel=col_vec([0,10,0]*1e-3); 
sim_opt.osc.amp=col_vec(osc_amp_vel./sim_opt.osc.omega); 
sim_opt.osc.lambda=0;% 0.30;%5/500;
sim_opt.osc.phase=col_vec([0,0,(1/4)])*pi;
sim_opt.grav=const.g0;
sim_opt.grav=const.g0;
sim_opt.nsamp=3.2e4;%number of atoms to simulate
sim_opt.verbose=0;

times=linspace(0,3*(2*pi/sim_opt.osc.omega(2)),50);
sim_opt.progress_bar=false;
jjmax=numel(times);


res_st={};
%res_st.atom_num=logspace(log10(500),log10(1e8),20)';
res_st.atom_num=logspace(log10(10),log10(500),20)';
res_st.atom_num=res_st.atom_num(randperm(numel(res_st.atom_num)))
res_st.phase_shift.val=nan*res_st.atom_num;
res_st.phase_shift.se=nan*res_st.atom_num;
iimax=numel(res_st.atom_num);
for ii=1:iimax
    tf_param.num=res_st.atom_num(ii);

    vel_compare=[];
    vel_compare.bec.mean=nan(jjmax,3);
    vel_compare.pal.mean=nan(jjmax,3);
    vel_compare.pal.std=nan(jjmax,3);
    vel_compare.pal.ste=nan(jjmax,3);
    vel_compare.time=times;
    
    sim_outs=cell(jjmax,1);
    parfor jj=1:jjmax
        sim_opt_tmp=sim_opt;
        % do the simulation of a out-coupling time
        fprintf('\n a. num %u of %u, time %u of %u \n',ii,iimax,jj,jjmax) 
        tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);
        
        sim_opt_tmp.osc.time_offset=times(jj);
        sim_out=class_sim_outcoupling(tf_details,sim_opt_tmp);
        sim_outs{jj}=sim_out;
    end
    for jj=1:jjmax
        sim_out=sim_outs{jj};
        vel_compare.bec.mean(jj,:)=row_vec(sim_out.bec.vel);
        vel_compare.pal.mean(jj,:)=sim_out.final.vel.mean;
        vel_compare.pal.std(jj,:)=sim_out.final.vel.std;
        vel_compare.pal.ste(jj,:)=sim_out.final.vel.std/...
                                    sqrt(size(sim_out.final.vel.vals,1));
    end
    
    fit_idx=2;
    data_in=[];
    data_in.x=vel_compare.pal.mean(:,fit_idx);
    data_in.t=vel_compare.time;
    data_in.unc_x=vel_compare.pal.ste(:,fit_idx);
    opts_in=[];
    fit_pal=fit_sine_to_data(data_in,opts_in);
    data_in=[];
    data_in.x=vel_compare.bec.mean(:,fit_idx);
    data_in.t=vel_compare.time;
    fit_bec=fit_sine_to_data(data_in,opts_in);
    % phase is in multiples of 2 pi
    res_st.phase_shift.val(ii)=fit_pal.fitparam{'phase','Estimate'}-fit_bec.fitparam{'phase','Estimate'};
    res_st.phase_shift.se(ii)=sqrt( fit_pal.fitparam{'phase','SE'}^2+fit_bec.fitparam{'phase','SE'}^2);

end

%%
%load('obs_phase_shift_with_anum.mat','res_st')
font_name='cmr10';
linewidth=1.5;
font_size=15;

stfig('velocity error')
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

colors_main=[[210,56,46];[50,163,73];[73,133,173]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2)-0.2;
colors_shaded=colorspace('HSV->RGB',hsv);


stfig('phase shift')
yval=rad2deg(wrapToPi(res_st.phase_shift.val*2*pi));
yval(isoutlier(yval,'mean'))=nan;
yerr=rad2deg(res_st.phase_shift.se*2*pi);
color_idx=3;
errorbar(res_st.atom_num,yval,yerr,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(color_idx,:),...
     'LineWidth',1.5,'MarkerFaceColor',colors_shaded(color_idx,:))
y_mean=nanmean(yval);
y_std=nanstd(yval);
yline(y_mean,'Color',colors_shaded(1,:),'LineWidth',linewidth)
yline(y_mean+y_std,'--','Color',colors_shaded(2,:),'LineWidth',linewidth)
yline(y_mean-y_std,'--','Color',colors_shaded(2,:),'LineWidth',linewidth)
legend('simulation','mean','$\pm\sigma$ ')
ln=legend('FontSize',font_size*0.9);
legend('location','best')
ln.EdgeColor=[1,1,1];
xlabel('Atom Number')
ylabel('Phase Shift (deg)')
set(gca, 'XScale', 'log')
xlim([100,2e8])
ylim([-5.6,-2.5])
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
set(gcf,'Position',[-1800,355,676,453])
set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])
%grid on
fig_name='pal_apparent_phase_shift';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
exportgraphics(gca,fullfile(fig_dir,strcat(fig_name,'.pdf')))
