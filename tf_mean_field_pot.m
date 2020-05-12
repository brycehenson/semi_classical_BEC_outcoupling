function [potential,pgrad,p2grad]=tf_mean_field_pot(pos_sample,tf_param,attenuation)

if attenuation~=1
    error('not yet implmeneted')
end
if nargout==1
    trap_potential=pot_harmonic_trap(pos_sample,tf_param);
    potential=(tf_param.mu_chem_pot-trap_potential);
    potential(potential<0)=0;
elseif nargout==2
    [trap_potential,pgrad]=pot_harmonic_trap(pos_sample,tf_param);
    potential=(tf_param.mu_chem_pot-trap_potential);
    out_bec_mask=potential<0;
    potential(out_bec_mask)=0;
    pgrad=-pgrad;
    pgrad(out_bec_mask,:)=0;
elseif nargout==3
    [trap_potential,pgrad,p2grad]=pot_harmonic_trap(pos_sample,tf_param);
    potential=(tf_param.mu_chem_pot-trap_potential);
        out_bec_mask=potential<0;
    potential(out_bec_mask)=0;
    pgrad=-pgrad;
    pgrad(out_bec_mask,:)=0;
    p2grad=-p2grad;
    p2grad(out_bec_mask,:)=0;
end





end