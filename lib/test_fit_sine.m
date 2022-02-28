%test_fit_sine

t_samp=linspace(0,10,1e4);

d_sine=[];
d_sine.omega=10;
x_val=deriv_damp_sine_wave(t_samp,d_sine);
x_err=x_val*0+0.01;

plot(t_samp,x_val)

data_in=[];
data_in.x=x_val;
data_in.t=t_samp;
opts_in=[];
fit_res=fit_sine_to_data(data_in,opts_in)



%%
fit_osc_plot_handle='test1'
font_name='cmr10';
font_size_global=10;
font_size_label=10;
colors_main=[[255,0,0];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);


tplotvalues=linspace(min(t_samp),max(t_samp),1e5)';
[prediction,ci]=predict(fit_res.fitobject,tplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
stfig(fit_osc_plot_handle);
clf;


subplot(2,1,1)
plot(t_samp,x_val,'kx-')
ylabel('V_{x} (mm/s)')
xlabel('Time (s)')
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
legend('x','y','z')

subplot(2,1,2)
time_start=min(t_samp);

shaded_ci_lines=false;

if shaded_ci_lines
    patch([predictorplot(:,1)', fliplr(predictorplot(:,1)')]-time_start, [ci(:,1)', fliplr(ci(:,2)')]*1e3, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
    hold on
else
    plot(tplotvalues-time_start,ci(:,1),'-','LineWidth',1.5,'Color',color_shaded)
    hold on
    plot(tplotvalues-time_start,ci(:,2),'-','LineWidth',1.5,'Color',color_shaded)
end  
plot(tplotvalues-time_start,prediction,'-','LineWidth',1.0,'Color',colors_main(3,:))
ax = gca;
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
errorbar(t_samp-time_start,x_val,x_err,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5) 
set(gcf,'Color',[1 1 1]);
xlabel('Time (s)','FontSize',font_size_label)
ylabel('V_{x} (mm/s)','FontSize',font_size_label)
title(sprintf('amp=%.2f±%.2f mm/s,omega=%.2f±%.2f Hz,Damp=%.2f±%.2f s',...
     fit_res.fitobject.Coefficients.Estimate(1)*1e3, fit_res.fitobject.Coefficients.SE(1)*1e3,...
      fit_res.fitobject.Coefficients.Estimate(3), fit_res.fitobject.Coefficients.SE(3),...
      fit_res.fitobject.Coefficients.Estimate(2), fit_res.fitobject.Coefficients.SE(2)))
hold off
ax = gca;
xlim([tplotvalues(1)-0.01,tplotvalues(end)+0.01]-time_start)
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth',1.0)
set(gca,'FontSize',font_size_global,'FontName',font_name)
%saveas(gca,sprintf('%sfit_dld_shot_num%04u.png',anal_opts_osc_fit.global.out_dir,dld_shot_num))




%% 


font_name='cmr10';
font_size_global=10;
font_size_label=10;
colors_main=[[255,0,0];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);


tplotvalues=linspace(min(txyz_tmp(:,1)),...
    max(txyz_tmp(:,1)),1e5)';
predictorplot=[tplotvalues,...
           interp1(predictor(:,1),predictor(:,2),tplotvalues),...
           interp1(predictor(:,1),predictor(:,3),tplotvalues)];
[prediction,ci]=predict(fitobject,predictorplot,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
stfig(fit_osc_plot_handle);
clf;


subplot(2,1,1)
plot(txyz_tmp(:,1),txyz_tmp(:,2)*1e3,'kx-')
hold on
plot(txyz_tmp(:,1),txyz_tmp(:,3)*1e3,'rx-')
plot(txyz_tmp(:,1),txyz_tmp(:,4)*1e3,'bx-')
hold off
ylabel('V_{x} (mm/s)')
xlabel('Time (s)')
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
legend('x','y','z')

subplot(2,1,2)
time_start=min(predictor(:,1));

shaded_ci_lines=false;

if shaded_ci_lines
    patch([predictorplot(:,1)', fliplr(predictorplot(:,1)')]-time_start, [ci(:,1)', fliplr(ci(:,2)')]*1e3, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
    hold on
else
    plot(predictorplot(:,1)-time_start,ci(:,1)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
    hold on
    plot(predictorplot(:,1)-time_start,ci(:,2)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
end  
plot(predictorplot(:,1)-time_start,prediction*1e3,'-','LineWidth',1.0,'Color',colors_main(3,:))
ax = gca;
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
errorbar(predictor(:,1)-time_start,txyz_tmp(:,anal_opts_osc_fit.dimension+1)*1e3,xyzerr_tmp(:,anal_opts_osc_fit.dimension)*1e3,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5) 
set(gcf,'Color',[1 1 1]);
xlabel('Time (s)','FontSize',font_size_label)
ylabel('V_{x} (mm/s)','FontSize',font_size_label)
title(sprintf('amp=%.2f±%.2f mm/s,omega=%.2f±%.2f Hz,Damp=%.2f±%.2f s',...
    fitobject.Coefficients.Estimate(1)*1e3,fitobject.Coefficients.SE(1)*1e3,...
     fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2),...
     fitobject.Coefficients.Estimate(7),fitobject.Coefficients.SE(7)))
hold off
ax = gca;
xlim([predictorplot(1,1)-0.01,predictorplot(end,1)+0.01]-time_start)
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth',1.0)
set(gca,'FontSize',font_size_global,'FontName',font_name)
saveas(gca,sprintf('%sfit_dld_shot_num%04u.png',anal_opts_osc_fit.global.out_dir,dld_shot_num))


%% do a speprate plot of the fit for the paper     
stfig(fit_osc_nice_plot_handle)
clf
shaded_ci_lines=false;

if shaded_ci_lines
    patch([predictorplot(:,1)', fliplr(predictorplot(:,1)')]-time_start, [ci(:,1)', fliplr(ci(:,2)')]*1e3, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
    hold on
else
    plot(predictorplot(:,1)-time_start,ci(:,1)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
    hold on
    plot(predictorplot(:,1)-time_start,ci(:,2)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
end  
plot(predictorplot(:,1)-time_start,prediction*1e3,'-','LineWidth',1.0,'Color',colors_main(3,:))
prev_gca_pos=get(gca,'Position');
set(gca,'Position',[0.1,0.2,0.88,0.78])
ax = gca;
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
errorbar(predictor(:,1)-time_start,txyz_tmp(:,anal_opts_osc_fit.dimension+1)*1e3,xyzerr_tmp(:,anal_opts_osc_fit.dimension)*1e3,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5) 
set(gcf,'Color',[1 1 1]);
xlabel('Time (s)','FontSize',font_size_label)
yticks(-10:10:10)
ylim([-16,16])
ylabel('$V_{x}$ (mm/s)','FontSize',font_size_label)
label_text=sprintf('amp=%.2f±%.2f mm/s,omega=%.2f±%.2f Hz,Damp=%.2f±%.2f s',...
    fitobject.Coefficients.Estimate(1)*1e3,fitobject.Coefficients.SE(1)*1e3,...
     fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2),...
     fitobject.Coefficients.Estimate(7),fitobject.Coefficients.SE(7));
an=annotation('textbox',[.3 .3 .5 .6],'String',label_text,'FitBoxToText','on');
an.FontSize = 5;
hold off
ax = gca;
xlim([predictorplot(1,1)-0.01,predictorplot(end,1)+0.01]-time_start)
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth',1.0)
set(gca,'FontSize',font_size_global,'FontName',font_name)
set(gcf,'Position',[100,100,270*3,100*3])
saveas(gca,sprintf('%sfit_nice_dld_shot_num%04u.png',anal_opts_osc_fit.global.out_dir,dld_shot_num))
export_fig(sprintf('%sfit_nice_dld_shot_num%04u.svg',anal_opts_osc_fit.global.out_dir,dld_shot_num))
export_fig(sprintf('%sfit_nice_dld_shot_num%04u.eps',anal_opts_osc_fit.global.out_dir,dld_shot_num))
%%
stfig(fit_resid_plot_handle);
subplot(4,1,1)
[fit_model_vals,fit_model_ci]=predict(fitobject,predictor,'Alpha',1-erf(1/sqrt(2)));
fit_model_ci=range(fit_model_ci,2)./2;
fit_resid=response-fit_model_vals;
fit_resid_sig=fit_resid./(sqrt(err_fit_dim.^2+fit_model_ci.^2));
plot(predictor(:,1),fit_resid,'k')
ylabel('residuals')
subplot(4,1,2)
plot(predictor(:,1),fit_resid_sig,'k')
ylabel('residuals sigma')
subplot(4,1,3)
fft_resid=fft_tx(predictor(:,1),fit_resid,'padding',10);
plot(fft_resid(1,:),abs(fft_resid(2,:)),'k')
hold on
plot(fft_resid(1,:),imag(fft_resid(2,:)),'b')
plot(fft_resid(1,:),real(fft_resid(2,:)),'r')
hold off
ylabel('residual spectral amplitude')
