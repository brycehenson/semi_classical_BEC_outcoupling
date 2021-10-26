function [pos,dpos]=osc_trap_pos(t,osc_det)

pos=damp_sine_wave(t,osc_det);
dpos=deriv_damp_sine_wave(t,osc_det);


end