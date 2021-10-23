function [pos,dpos]=osc_trap_pos(t,osc_det)

tau=osc_det.decay_tau;
amp=osc_det.amp;
omega=osc_det.omega;
phi=osc_det.phase;
time_offset=osc_det.time_offset;

t=t+time_offset;


lambda=1/tau;

damp_ratio=lambda./sqrt((lambda^2)+(omega));
omega_dampened=sqrt(1-(damp_ratio.^2)).*omega;


pos=amp.*exp(-t/tau).*sin(t*omega_dampened+phi);

% calculate the derivative
%if nargout
% beacuase no a pure sine wave we get a small phase shift compared to
% the pi/2 expected for the undamped case
%damp_phase_shift=asin(sqrt(1-(damp_ratio^2)));
% https://math.stackexchange.com/questions/2849945/sum-of-sin-waves-with-same-frequency-and-different-amplitudes-and-phase
A1=omega_dampened;
phi1=pi/2+phi;
A2=-lambda;
phi2=phi;
damp_phase_shift=atan( ( A1.*sin(phi1)+A2.*sin(phi2) ) ./...
                        ( A1.*cos(phi1)+A2.*cos(phi2) ) );
amp_fac=sqrt(omega_dampened.^2  + lambda.^2 );
dpos=-amp.*exp(-t/tau).*amp_fac.*sin(t*omega_dampened+phi+damp_phase_shift );
%end

end