function out_st=fit_sine_to_data(data,opts)

% INPUTS
% data.t 
% data.x
% data.unc_x

% opts.guess.freq

% opts.global_fit.enable
% opts.global_fit.param_range

%%

%warning('off','stats:nlinfit:ModelConstantWRTParam');
%warning('off','MATLAB:rankDeficientMatrix');
out_st=[]; %clear the output struct


if ~isfield(data,'unc_x')
    data.unc_x=[];
end
if ~isfield(opts,'guess')
    opts.guess=struct();
end

if ~isfield(opts.guess,'amp')
    opts.guess.amp='adaptive';
end
if ~isfield(opts.guess,'damp')
    opts.guess.damp=0;
end
if ~isfield(opts.guess,'freq')
    opts.guess.freq='adaptive';
end
if ~isfield(opts.guess,'phase')
    opts.guess.phase='adaptive';
end
if ~isfield(opts.guess,'offset')
    opts.guess.offset='adaptive';
end
if ~isfield(opts.guess,'grad')
    opts.guess.grad='adaptive';
end


        
%%
%try to find the peak osc freq to start the fit there
if strcmp(opts.guess.freq,'adaptive') ||  strcmp(opts.guess.phase,'adaptive') || strcmp(opts.guess.amp,'adaptive')
    dom_opt=[];
    dom_opt.num_components=1;
    dom_out=dominant_freq_components(data.t ,data.x,dom_opt);
    if strcmp(opts.guess.freq,'adaptive')
        opts.guess.freq=dom_out.freq(1);
    end
    if strcmp(opts.guess.phase,'adaptive')
        opts.guess.phase=dom_out.phase(1);
    end
    if  strcmp(opts.guess.amp,'adaptive')
         opts.guess.amp=dom_out.amp(1)*1.2;
    end
end

if strcmp(opts.guess.offset,'adaptive')
    opts.guess.offset=mean(data.x);
end

if strcmp(opts.guess.grad,'adaptive')
    len_dat=size(data.t,2);
    start_range=[1,round(len_dat/3)];
    end_range=[len_dat-round(len_dat/3):len_dat];
    grad_est=(mean(data.x(end_range(1):end_range(2)))-mean(data.x(start_range(1):start_range(2))))/...
        (data.t(end_range(1))-data.t(start_range(1)));
    opts.guess.grad=grad_est;
end


       
        
%%
predictor=data.t;
unc_response=data.unc_x;
if isempty(unc_response)
    unc_response=ones(size(predictor));
end
weights=1./(unc_response.^2);
weights(isnan(weights))=1e-20; %if nan then set to smallest value you can
weights=weights/sum(weights);
%
response=data.x;
    % the model is a standard dampened sine wave
    % add a max function to the damping rate so that if the fit model sets it
    % to a rediculous -ve number then we dont get an inf
    modelfun = @(b,x) b(1).*exp(-x.*max(-5e2,b(2))).*sin(b(3)*x*pi*2+b(4)*pi*2)+b(5)+b(6)*x;


if  (isfield(opts,'global_fit')) && opts.global_fit.enable
    %% do a global fit to get a good place to start the fit routine
    % before we call fitnlm we will try to get close to the desired fit parameters using a robust global
    % optimizer
    % simple model
    %modelfun_simple = @(b,x) exp(-x(:,1).*max(0,b(6))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(5)*x(:,1);
    gf_opt=[];
    gf_opt.domain=[opts.global_fit.param_range.amp;...  
                   opts.global_fit.param_range.damp;...     
                   opts.global_fit.param_range.freq;...     
                   opts.global_fit.param_range.phase;... 
                   opts.global_fit.param_range.offset;...   
                   opts.global_fit.param_range.grad;... 
                   ];        
    gf_opt.start=[fit_amp, fit_freq, fit_phase, fit_offset,0,0,2,0];
    gf_opt.rmse_thresh=2e-3;
    gf_opt.plot=false;
    gf_out=global_fit(predictor,response,modelfun,gf_opt);
else
    gf_out.params=[opts.guess.amp, opts.guess.damp, opts.guess.freq, opts.guess.phase,opts.guess.offset,opts.guess.grad];
end
%% set up the complex model with xterms


cof_names={'amp','damp','freq','phase','offset','grad'};
opt = statset('TolFun',1e-10,'TolX',1e-10,...
    'MaxIter',1e4,... %1e4
    'UseParallel',1);
beta0=gf_out.params;
         
%%
fitobject=fitnlm(predictor,response,modelfun,beta0,...
    'Weights',weights,'options',opt,...
    'CoefficientNames',cof_names);
out_st.fitobject=fitobject;
out_st.fitparam=fitobject.Coefficients;
out_st.model_coefs=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];
out_st.fit_rmse=fitobject.RMSE;


% we will later use this fit to find the undamped frequency of the oscillator
%limiting frequnecy prediction from http://adsabs.harvard.edu/full/1999DSSN...13...28M
% meanwidth=sqrt(mean(squeeze(data.mcp_tdc.al_pulses.vel_zxy.std(ii,:,anal_opts_osc_fit.dimension)).^2))*1e3;
% frequnclim=sqrt(6/sum(data.mcp_tdc.al_pulses.num_counts(ii,:)))*...
%     (1/(pi*range(data.t)))*...
%     (meanwidth/fitparam{2,1});
%fprintf('sampling limit %2.3g Hz, fint unc %2.3g Hz, ratio %2.3g \n',[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim])
% out_st.fit_sample_limit{ii}=[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim];


        
end
    