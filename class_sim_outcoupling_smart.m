function out_st=class_sim_outcoupling_smart(tf_details,opt)
% this is a dump of some code which was an attempt to make a smarter outcoupling simulation
% ultimately i dont think its that worth it, at least how this code is written
% i think a better approach is for the stationary bec case is to solve 
% for when the cosh(t omega) trajectory in each axis reaches the TF radius
% then convert to velocity
% and only sim those atoms which are going up at an angle that will cobe back to intersect
% but im not super interested in this case at the moment


chunk_size=1000;
data_in=pos_sample;
tprop_all=5e-3;
%find the time after which all the particles are in freefall
% first we find the maximum velocity upwards by equating the chem potential to the kenetic energy up
% (1/2)·m·vmax^2=tf_details.mu_chem_pot;
vmax =sqrt(tf_details.mu_chem_pot*2/tf_details.inputs.mass);
%time for velocity to go to zero is vmax/g, and to come back to BEC is twice this
% a simple formula would be
% tprop_end=2*vmax/g;
% given y(0)=R_tf find y(t)=-R_tf
% -R_tf =vmax -(1/2)g t^2 +R_tf
tprop_end=sqrt(2*(2*tf_details.tf_radi(3) +vmax)/const.g0);
tprop_end=tprop_end*0.25; %hack for testing


time_sim=tic;
outcoupled.pos=nan*pos_sample;
outcoupled.vel=nan*pos_sample;
iimax=ceil(size(data_in,1)/chunk_size);
iimax=10;
ode_opts = odeset('InitialStep',1e-4,'MaxStep',1e-3);
fprintf('\nrunning simulation\n')
output_chars=fprintf('chunk %04u of %04u',0,iimax);

for ii=1:iimax
    start_idx=(ii-1)*chunk_size+1;
    end_idx=start_idx+chunk_size-1;
    end_idx=min([end_idx,size(data_in,1)]);
data_chunk=data_in(start_idx:end_idx,:);
    xstart=nan([size(data_chunk),2]);
    xstart(:,:,1)=data_chunk;
    xstart(:,:,2)=0;
    xshape=size(xstart);

    
    xsol=ode45(@(t,x) state_update_fun_grav(x,tf_details,xshape,const.g0),[0,tprop_all],...
        reshape(xstart,[],1),...
        ode_opts);
    % lets find the position and velocity at the final time
    xoutc=deval(xsol,tprop_all);
    xoutc=reshape(xoutc,xshape);
    
    %check that no atoms are still in the condensate (or rather a cube approx therof)
    in_condensate_mask=abs(xoutc(:,:,1))<repmat(tf_details.tf_radi',[size(xoutc,1),1]);
    in_condensate_mask=all(in_condensate_mask,2);
    num_in_condensate=sum(in_condensate_mask);
    if num_in_condensate>0
        fprintf('\n%u atoms still in condensate \n', num_in_condensate)
        fprintf(repmat(' ',[1,output_chars]))
    end
    
    % lets mask out the particles that require further simulation
    % to start with we just select everything thats above z=0
    second_sim_mask=xoutc(:,3,1) >0;
    % find the time that the atoms fall through z=0
    % z(t)=z0+t*v_z+(1/2) a_z t^2
    %t_crossing= (vz + sqrt(vz^2+2*g*y0) )/g
    t_crossing=( xoutc(second_sim_mask,3,2) +...
                 sqrt( xoutc(second_sim_mask,3,2).^2+2*const.g0*xoutc(second_sim_mask,3,1) )...
                )./const.g0;
    pos_crossing=xoutc(second_sim_mask,1:2,1)+xoutc(second_sim_mask,1:2,2).*repmat(t_crossing,[1,2]);
    % lets see if x&y are within the tf radius at the crossing
    second_sim_mask(second_sim_mask)=all(abs(pos_crossing)<repmat(tf_details.tf_radi(1:2)',[size(pos_crossing,1),1]),2);
    
    % those particles that are not going up we propagate to t=tprop_end
    tdelta=tprop_end-tprop_all;
    % we use basic kinematics x(t)=x0+t*v+(1/2) a t^2
    xfinal=xoutc;
    xfinal(second_sim_mask,:,:)=nan;
    % increment the position using the current velocity
    xfinal(~second_sim_mask,:,1)=xfinal(~second_sim_mask,:,1)+xfinal(~second_sim_mask,:,2)*tdelta;
    % add the (1/2) a t^2 component in z
    xfinal(~second_sim_mask,3,1)=xfinal(~second_sim_mask,3,1)-(1/2)*const.g0*(tdelta^2);
    % increment the velocity in z
    xfinal(~second_sim_mask,3,2)=xfinal(~second_sim_mask,3,2)-const.g0*tdelta;
    
    if sum(second_sim_mask)>0
        x_more_sim=xoutc(second_sim_mask,:,:);
        xshape=size(x_more_sim);
        xsol=ode45(@(t,x) state_update_fun_grav(x,tf_details,xshape,const.g0),[0,tdelta],...
            reshape(x_more_sim,[],1),...
            ode_opts);
        % lets find the position and velocity at the final time
        xout_second=deval(xsol,tdelta);
        xout_second=reshape(xout_second,xshape);
        xfinal(second_sim_mask,:,:)=xout_second;
    end
    
    
    outcoupled.pos(start_idx:end_idx,:)=xfinal(:,:,1);
    outcoupled.vel(start_idx:end_idx,:)=xfinal(:,:,2);
    fprintf(repmat('\b',[1,output_chars]))
    output_chars=fprintf('chunk %04u of %04u',ii,iimax);
end

toc(time_sim) 



end