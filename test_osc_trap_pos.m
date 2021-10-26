% test_osc_trap_pos

do_plots=true;


inputs=[];
inputs.lambda=1/5;
inputs.amp=8.5726;
inputs.omega=12.2846;
inputs.phase=-0.6*pi;%-pi*1.1;
inputs.time_offset=0;

t_vec=linspace(0,20,1e5);
t_vec=col_vec(t_vec);
[pos_vec,deriv_vec]=osc_trap_pos(t_vec,inputs);
fun_pointer=@(t) osc_trap_pos(t,inputs);
num_deriv=num_grad_vecs(fun_pointer,t_vec,1e-6,true);

if do_plots
    stfig('testing_osc')
    clf
    subplot(3,1,1)
    plot(t_vec,pos_vec)
    subplot(3,1,2)
    plot(t_vec,deriv_vec,'k')
    %num_deriv=derivest(fun_pointer,t_vec);
    hold on
    plot(t_vec,num_deriv,'r')
    hold off
end

deriv_err=deriv_vec-num_deriv;

if do_plots
    subplot(3,1,3)
    plot(t_vec,deriv_err,'k')
end

rms_deriv_err=rms(deriv_err);
rms_amp=rms(pos_vec);

% if rms_deriv_err>1e-9
%     error('too large')
% end
fprintf('rms err %G\n',rms_deriv_err)
fprintf('norm rms err %G\n',rms_deriv_err/(rms_amp*inputs.omega))

%%
for ii=1:100
   test_osc_trap_pos() 
end


function test_osc_trap_pos()
    inputs=[];
    inputs.lambda=1/rand_interval([1,100],1);
    inputs.amp=rand_interval([1,10],1);
    inputs.omega=rand_interval([1,100],1);
    inputs.phase=rand_interval([-5,5],1)*pi;
    inputs.time_offset=rand_interval([-100,100],1);

    t_vec=linspace(0,20,1e5);
    t_vec=col_vec(t_vec);
    [pos_vec,deriv_vec]=osc_trap_pos(t_vec,inputs);
    fun_pointer=@(t) osc_trap_pos(t,inputs);
    num_deriv=num_grad_vecs(fun_pointer,t_vec,1e-6*(1/inputs.omega),true);


    deriv_err=deriv_vec-num_deriv;
    rms_amp=rms(pos_vec);
    rms_deriv_err=rms(deriv_err)/(rms_amp*inputs.omega);
    if rms_deriv_err>1e-6
        rms_deriv_err
        inputs
        error('error too large %g',rms_deriv_err)
    end
end