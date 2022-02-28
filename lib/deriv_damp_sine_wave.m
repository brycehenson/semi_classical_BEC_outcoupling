function dpos=deriv_damp_sine_wave(t,params)
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
% calculate the derivative
% beacuase not a pure sine wave we get a small phase shift compared to
% the pi/2 expected for the undamped case
%damp_phase_shift=asin(sqrt(1-(damp_ratio^2)));
% https://math.stackexchange.com/questions/2849945/sum-of-sin-waves-with-same-frequency-and-different-amplitudes-and-phase
% https://dspguru.com/files/Sum_of_Two_Sinusoids.pdf
% https://lpsa.swarthmore.edu/BackGround/phasor/phasor.html
% We can use a trig identity (a·cos(x)+b·sin(x)) = c·cos(x+φ) where c=√a2+b2 and φ=-atan2(b,a).
% https://math.stackexchange.com/questions/2085687/proving-a-cos-xb-sin-x-sqrta2b2-cosx-alpha

A1=omega_dampened;
phi1=phi;
A2=-lambda;
phi2=phi;

damp_phase_shift=-atan( -lambda./omega_dampened);
amp_fac=sqrt(omega_dampened.^2  + lambda.^2 );
dpos=amp.*exp(-t*lambda).*amp_fac.*cos(t*omega_dampened+phi+damp_phase_shift );
%end

end