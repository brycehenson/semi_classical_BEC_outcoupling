function x_val=damp_sine_wave(t,params)
% params.decay_lambda
% params.amp
% params.omega
% params.phase
% params.time_offset


if ~isfield(params,'lambda')
    params.lambda=0;
end
if ~isfield(params,'time_offset')
    params.time_offset=0;
end
if ~isfield(params,'phase')
    params.phase=0;
end
if ~isfield(params,'amp')
    params.amp=1;
end
if ~isfield(params,'omega')
    params.omega=1;
end

lambda=params.lambda;
amp=params.amp;
omega=params.omega;
phi=wrapToPi(params.phase);
time_offset=params.time_offset;

t=t+time_offset;
damp_ratio=lambda./sqrt((lambda^2)+(omega));
omega_dampened=sqrt(1-(damp_ratio.^2)).*omega;
x_val=amp.*exp(-t.*lambda).*sin(t*omega_dampened+phi);

end