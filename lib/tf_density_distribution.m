function density=tf_density_distribution(pos,tf_param)

potential=pot_harmonic_trap(pos,tf_param);
density=(tf_param.mu_chem_pot-potential)./tf_param.u_eff_interaction;
density(density<0)=0;


end