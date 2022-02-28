%test_tf_mean_field_pot
tf_p_in=[];
tf_p_in.trap_freq=2*pi*[55,411,415]; 
tf_p_in.mass=const.mhe;
tf_p_in.a_scat_len=const.ahe_scat;
tf_p_in.num=3e5;
tf_details=bec_properties(tf_p_in.trap_freq/(2*pi),tf_p_in.num,tf_p_in.mass,tf_p_in.a_scat_len);


nsamp=10;
% sample in a box 1.05x the tf radi
pos_sample=2*(rand(nsamp,3)-0.5);
pos_sample=pos_sample.*repmat(tf_details.tf_radi',[10,1])*1.05;

[trap_potential,pgrad]=tf_mean_field_pot(pos_sample,tf_details,1);

delt=tf_details.tf_radi_bar*1e-6;
ngrad=num_grad_vecs(@(x) tf_mean_field_pot(x,tf_details,1),pos_sample,delt,1);

graderr=ngrad-pgrad;
mean_frac_graderr=mean(graderr(:))/mean(pgrad(:))

if abs(mean_frac_graderr)>1e-6
    error('numeric derivative is different to analytical')
end
