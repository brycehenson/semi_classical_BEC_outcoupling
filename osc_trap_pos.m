function [pos,dpos]=osc_trap_pos(t,osc_det)

pos=osc_det.amp.*sin(t*osc_det.omega+osc_det.phase);

if nargout
    dpos=osc_det.omega.*osc_det.amp.*cos(t*osc_det.omega+osc_det.phase);
end

end