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

if ~isfield(opt,'den_power')
    opt.den_power=1;
end




% (1/2) m v^2= 1/2 m omega^2 x^2
if opt.verbose>1
    det_time=0.4;
    fprintf('osc det displacement %s mm\n',sprintf('%.2f,',opt.osc.amp.*opt.osc.omega*det_time*1e3))
end



if isfield(opt,'t_end')
    tprop_end=opt.t_end; 
else
    % find the outcoupling time
    % cosh(t omega)*r(0)=rtf
    r_start=0.01; % radius at which the atom starts
    % the time to exit can be given as
    t_out_couple=acosh(1/r_start)/max(opt.osc.omega); 
    tprop_end=t_out_couple;
    if opt.grav>0
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
    end
end
fprintf('simulation t_end=%.2g \n',tprop_end)

[bec_pos,bec_vel]=osc_trap_pos(0,opt.osc);

expanding_bec=false;
if isfield(opt,'trap_off') && ~isnan(opt.trap_off) && opt.trap_off<tprop_end
    expanding_bec=true;
    % find the scaling solution function
    [~,~,lambda_fun]=tf_expand_scaling_trap_off_num(tf_details.inputs.omega,tprop_end-opt.trap_off);
    lambda_fun_shifted=@(t) lambda_fun(t-opt.trap_off);
    fprintf('lambda at end %s \n',sprintf('%.2f,',lambda_fun_shifted(tprop_end)))
else
    lambda_fun_shifted=@(t) nan; %hack to get parfor to work
end
% for debug
% plot(inlidx(lambda_fun(linspace(0,0.00001,1e3)),{':',1}))


time_sim=tic;

% adaptive chunk size
% too low a chunk size will slow things down from more communication
% overhead, but if you ar doing sims with low numbers it can mean that not
% all the threads are used
min_chunk_size=100;
max_chunk_size=10000;
chunk_size=bound(round(opt.nsamp/get_num_workers()),min_chunk_size,max_chunk_size);
%chunk_size=100;

iimax=round(opt.nsamp/chunk_size);
fprintf('\nrunning simulation\n')
ode_opts = odeset('RelTol',1e-4,...
                    'InitialStep',0.01/max(tf_details.inputs.omega),...
                    'MaxStep',0.1/max(tf_details.inputs.omega),...
                    'Vectorized','off');
                    % the jacobian does not seem to speed things up
                    % 'Jacobian',@(t,x) state_update_jacobian(x,xshape,tf_details)

tmp_fin_states=cell(1,iimax);
tmp_start_states=cell(1,iimax);

% set up the parfor progress bar
fprintf('\n\n')
parfor_progress_imp(iimax);

parfor ii=1:iimax
    data_chunk=sample_pts_from_tf(chunk_size,tf_details,opt.den_power);
    % we define the state (pos,vel) of all these particles in a matrix
    % dimension = N (particle index) x 3(cart idx) x 2 (pos,vel)
    xstart=nan([size(data_chunk),2]);
    xstart(:,:,1)=data_chunk+repmat(bec_pos',[size(data_chunk,1),1]);
    % the velcoity will be the velocity of the bec at t=0
    xstart(:,:,2)=repmat(bec_vel',[size(data_chunk,1),1]);
    % optionally add the QD distribution
    if isfield(opt,'qd') 
        xstart(:,:,2)=xstart(:,:,2)+sample_from_k4_dist(chunk_size,opt.qd);
    end
    
    
    xshape=size(xstart);
    if expanding_bec
        xsol=ode45(@(t,x) state_update_fun_grav(x,tf_details,xshape,opt.grav,t,opt.osc,lambda_fun_shifted,opt.trap_off),[0,tprop_end],...
            reshape(xstart,[],1),...
            ode_opts);
    else
        %state_update_fun_grav(reshape(xstart,[],1),tf_details,xshape,opt.grav,1e-3,opt.osc)
        xsol=ode45(@(t,x) state_update_fun_grav(x,tf_details,xshape,opt.grav,t,opt.osc),[0,tprop_end],...
        reshape(xstart,[],1),...
        ode_opts);  
    end
    
    
    % lets find the position and velocity at the final time
    if xsol.x(end)==tprop_end %check that it found the final time
        xfinal=xsol.y(:,end);
    else
        xfinal=deval(xsol,tprop_end);
    end
    xfinal=reshape(xfinal,xshape);
    
    if expanding_bec
        lambda_end=lambda_fun_shifted(tprop_end);
    else
        lambda_end=[1,1,1];
    end
    %check that no atoms are still in the condensate (or rather a cube approx therof)
    lambda_end=repmat(lambda_end,[xshape(1),1]);
    xend_scaled=xfinal(:,:,1)./lambda_end;
    in_condensate_mask=abs(xend_scaled)<repmat(tf_details.tf_radi',[size(xfinal,1),1]);
    in_condensate_mask=all(in_condensate_mask,2);
    num_in_condensate=sum(in_condensate_mask);
    if num_in_condensate>0
        fprintf('\n%u atoms still in condensate \n', num_in_condensate)
        fprintf('\n')
        fprintf(repmat(' ',[1,50+12+36]))
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



function dxdt=state_update_fun_grav(xin,tf_param,xshape,grav,t,osc_det,lambda_fun,t_trap_off)
% this function returns the temporal derivatives of each of the state parameters in xin
% velocity and position for each cartesian axes are encoded in a vector
% for ease of processing the state is unfolded into a matrix
% the state is N (particle index) x 3(cart idx) x 2 (pos,vel)
% where xin(:,:,1) represents the postions of all particles
% where xin(:,:,2) represents the velocities
xin=reshape(xin,xshape);
dxdt=0*xin;
dxdt(:,:,1)=xin(:,:,2);

if nargin>6
    lambda_this=lambda_fun(t);
else
    lambda_this=[1,1,1];
    t_trap_off=inf;
end
lambds_prod=prod(lambda_this);
lambda_this=repmat(lambda_this,[xshape(1),1]);

if t<t_trap_off
    condensate_pos=osc_trap_pos(t,osc_det);
    condensate_pos(3)=condensate_pos(3)+(1/2)*grav*(t^2);
else
    % once the trap turns off the condensate will travel balisticaly
    [condensate_pos,dcondensate_pos]=osc_trap_pos(t_trap_off,osc_det);
    condensate_pos=condensate_pos+dcondensate_pos.*(t-t_trap_off);
    % if the trap has turned off
    condensate_pos(3)=condensate_pos(3)+(1/2)*grav*(t_trap_off^2)...
        - (1/2)*grav*((t-t_trap_off)^2) ;
end
pos_shifted=xin(:,:,1)-repmat(condensate_pos',[xshape(1),1]);
pos_scaled=pos_shifted./lambda_this;
[~,pot_grad]=tf_mean_field_pot(pos_scaled,tf_param,1);
% this scaling is from [Y. Castin and R. Dum, “Bose-Einstein Condensates in Time Dependent Traps,” Physical Review Letters, vol. 77, no. 27, pp. 5315–5319, Dec. 1996](https://doi.org/10.1103/PhysRevLett.77.5315) 
% the first scaling gives a derivative part of m \omega_j^2 r/\lambda_j
% now this scaling gives \omega_j^2 r/(\lambda_j^2 \prod(\lambda)) 
pot_grad_scaled=pot_grad./(lambda_this*lambds_prod);
force=-pot_grad_scaled;
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




