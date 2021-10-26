function unc_est=unc_est_in_pal_meas(tf_details,samp_param)

% inputs
% samp_param.pal_num_per_pulse
% samp_param.osc_amp
% samp_param.num_pulses
% samp_param.time

    %detected per pulse
    pal_num_per_pulse=samp_param.pal_num_per_pulse;
    num_pulses=samp_param.num_pulses;
    osc_amp=samp_param.osc_amp;
    time=samp_param.time;
    lambda=samp_param.damp_rate
    
    % estimate the uncertianty in a single pal pulse
    pos_unc_pal=sqrt(1/(12*pal_num_per_pulse))*2*tf_details.pal_vmax;
    
    in_st=[];
    in_st.sigma_obs=pos_unc_pal; %uncertianty in the observerd variable

    in_st.amp=osc_amp; % amplitude of the sine wave
    in_st.samp_num=num_pulses; % number of samples
    in_st.samp_time=time; %duration over which the oscillation was sampled
    in_st.damp_rate=lambda;
    
    unc_est=analy_err_in_fit_sine(in_st);
   

end