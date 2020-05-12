function out_st=class_sim_outcoupling(tf_details,opt)
% simulate the trajectory of atoms outcoupled from a TF bec
% the potential is entirely the mean field potential from the BEC
% no interations are considered in the atom laser
% this is a classical sim so interferences will not be treated properly
% inputs
%   tf_details
%   opts
%       .osc
%           .phase
%           .omega
%           .amp
%       .nsamp
%       .verbose

%TODO
% include start points in output structure

if ~isfield(opt,'verbose')
    opt.verbose=0;
end

chunk_size=200;


% (1/2) m v^2= 1/2 m omega^2 x^2
if opt.verbose>1
    det_time=0.4;
    fprintf('osc det displacement %s mm\n',sprintf('%.2f,',opt.osc.amp.*opt.osc.omega*det_time*1e3))
end

% find the outcoupling time
% cosh(t omega)*r(0)=rtf
r_start=0.01; % radius at which the atom starts
% the time to exit can be given as
t_out_couple=acosh(r_start)/max(opt.osc.omega); 
% note that there is always an unstable equlibrium where an atom stays in the condensate
% in the falling frame this is an atom 'surfing' on the bec potential and accelerating up at g

% find the time after which all the particles are in freefall
% first we find the maximum velocity upwards by equating the chem potential to the kenetic energy up
% (1/2)·m·vmax^2=tf_details.mu_chem_pot;
vmax =sqrt(tf_details.mu_chem_pot*2/tf_details.inputs.mass);
%time for velocity to go to zero is vmax/g, and to come back to BEC is twice this
% a simple formula would be
% tprop_end=2*vmax/g;
% given y(0)=R_tf find y(t)=-R_tf
% -2*R_tf =vmax t -(1/2)g t^2 
% we add on the oscillation amplitude here as well
tprop_end=( vmax+sqrt((vmax.^2)+4*(1/2)*opt.grav*2*(tf_details.tf_radi(3)+opt.osc.amp(3))) )/opt.grav;
% and add on a fudge factor
tprop_end=1.2*tprop_end;

[bec_pos,bec_vel]=osc_trap_pos(0,opt.osc);

time_sim=tic;

iimax=round(opt.nsamp/chunk_size);
fprintf('\nrunning simulation\n')
output_chars=fprintf('chunk %04u of %04u',0,iimax);
ode_opts = odeset('RelTol',1e-4,...
                    'InitialStep',0.01/max(tf_details.inputs.omega),...
                    'MaxStep',0.1/max(tf_details.inputs.omega),...
                    'Vectorized','off');
                    % the jacobian does not seem to speed things up
                    % 'Jacobian',@(t,x) state_update_jacobian(x,xshape,tf_details)

tmp_fin_states=cell(1,iimax);
tmp_start_states=cell(1,iimax);

% set up the parfor progress bar
parfor_progress_imp(iimax);

parfor ii=1:iimax
    data_chunk=sample_pts_from_tf(chunk_size,tf_details);
    % we define the state (pos,vel) of all these particles in a matrix
    % dimension = N (particle index) x 3(cart idx) x 2 (pos,vel)
    xstart=nan([size(data_chunk),2]);
    xstart(:,:,1)=data_chunk+repmat(bec_pos',[size(data_chunk,1),1]);
    xstart(:,:,2)=repmat(bec_vel',[size(data_chunk,1),1]);
    xshape=size(xstart);

    xsol=ode45(@(t,x) state_update_fun_grav(x,tf_details,xshape,opt.grav,t,opt.osc),[0,tprop_end],...
        reshape(xstart,[],1),...
        ode_opts);
    
    % lets find the position and velocity at the final time
    if xsol.x(end)==tprop_end %check that it found the final time
        xfinal=xsol.y(:,end);
    else
        xfinal=deval(xsol,tprop_end);
    end
    xfinal=reshape(xfinal,xshape);
    
    %check that no atoms are still in the condensate (or rather a cube approx therof)
    in_condensate_mask=abs(xfinal(:,:,1))<repmat(tf_details.tf_radi',[size(xfinal,1),1]);
    in_condensate_mask=all(in_condensate_mask,2);
    num_in_condensate=sum(in_condensate_mask);
    if num_in_condensate>0
        fprintf('\n%u atoms still in condensate \n', num_in_condensate)
        fprintf(repmat(' ',[1,output_chars]))
    end
    
    % we use this temp cell form so that the loop is parfor compatable 
    tmp_fin_states{ii}=xfinal;
    tmp_start_states{ii}=xstart;
%     fprintf(repmat('\b',[1,output_chars]))
%     output_chars=fprintf('chunk %04u of %04u',ii,iimax);
    parfor_progress_imp;
end
parfor_progress_imp(0);

start_pos=nan(iimax*chunk_size,3);
start_vel=nan(iimax*chunk_size,3);
outcoupled_pos=nan(iimax*chunk_size,3);
outcoupled_vel=nan(iimax*chunk_size,3);
for ii=1:iimax
    start_idx=(ii-1)*chunk_size+1;
    end_idx=start_idx+chunk_size-1;
    xstart=tmp_start_states{ii};
    start_pos(start_idx:end_idx,:)=xstart(:,:,1);
    start_vel(start_idx:end_idx,:)=xstart(:,:,2);
    xfinal=tmp_fin_states{ii};
    outcoupled_pos(start_idx:end_idx,:)=xfinal(:,:,1);
    outcoupled_vel(start_idx:end_idx,:)=xfinal(:,:,2);
    
end

time_to_sim=toc(time_sim);
time_for_single=time_to_sim/(iimax*chunk_size);

if opt.verbose>1
    fprintf('time for single particle sim %.2e \n',time_for_single)
end

out_st.start.pos.vals=start_pos;
out_st.start.pos.mean=mean(start_pos,1);
out_st.start.pos.std=std(start_pos,1);

out_st.start.vel.vals=start_vel;
out_st.start.vel.mean=mean(start_vel,1);
out_st.start.vel.std=std(start_vel,1);

out_st.final.pos.vals=outcoupled_pos;
out_st.final.pos.mean=mean(outcoupled_pos,1);
out_st.final.pos.std=std(outcoupled_pos,1);

out_st.final.vel.vals=outcoupled_vel;
out_st.final.vel.mean=mean(outcoupled_vel,1);
out_st.final.vel.std=std(outcoupled_vel,1);

out_st.final.time=tprop_end;
out_st.bec.vel=bec_vel;
out_st.bec.pos=bec_pos;
out_st.run_info.runtime=time_to_sim;
out_st.run_info.time_single=time_for_single;

end


% some test code for the below functions
% x_test=[1,1,1,0,0,0]
% x_test_shape=[1,3,2];
% x_deriv=state_update_fun(x_test,tf_details,x_test_shape)
% 
% %%
% data_chunk=pos_sample(1:4,:);
% xstart=nan([size(data_chunk),2]);
% xstart(:,:,1)=data_chunk;
% xstart(:,:,2)=0;
% xderiv_dirn=xstart*0;
% xshape=size(xstart);
% xstart=reshape(xstart,[],1)  ;
% 
% %[jac,jac_err] = jacobianest(@(x) state_update_fun_grav(x,tf_details,xshape,const.g0) ,xstart)
% fstart=state_update_fun_grav(xstart,tf_details,xshape,const.g0)



function dxdt=state_update_fun_grav(xin,tf_param,xshape,grav,t,osc_det)
% this function returns the temporal derivatives of each of the state parameters in xin
% velocity and position for each cartesian axes are encoded in a vector
% for ease of processing the state is unfolded into a matrix
% the state is N (particle index) x 3(cart idx) x 2 (pos,vel)
% where xin(:,:,1) represents the postions of all particles
% where xin(:,:,2) represents the velocities
xin=reshape(xin,xshape);
dxdt=0*xin;
dxdt(:,:,1)=xin(:,:,2);
condensate_pos=osc_trap_pos(t,osc_det);
condensate_pos(3)=condensate_pos(3)+(1/2)*grav*(t^2);
pos_shifted=xin(:,:,1)-repmat(condensate_pos',[size(xin,1),1]);
[~,pot_grad]=tf_mean_field_pot(pos_shifted,tf_param,1);
force=-pot_grad;
%f=ma;
accel=(force/tf_param.inputs.mass) ;
accel(:,3)=accel(:,3);
dxdt(:,:,2)=accel;
dxdt=reshape(dxdt,[],1);
end


function jac=state_update_jacobian(xin,xshape,tf_param)
% the jacobian representing how the state temporal derivative changes with the state variables
% is fairly sparse
nx=numel(xin);
jac=zeros([1,1]*numel(xin));

for ii=1:(nx/2)
    jac(ii,ii+(nx/2))=1;
end

xin=reshape(xin,xshape);
[~,~,pot_2grad]=tf_mean_field_pot(xin(:,:,1),tf_param,1);
for ii=1:(nx/2)
    jac(ii+(nx/2),ii)=-pot_2grad(ii);
end

end




