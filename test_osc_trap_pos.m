% test_osc_trap_pos

inputs=[];
inputs.decay_tau=5;
inputs.amp=1;
inputs.omega=1;
inputs.phase=0;

t_vec=linspace(0,20,1e5);
t_vec=col_vec(t_vec);
[pos_vec,deriv_vec]=osc_trap_pos(t_vec,inputs);
stfig('testing_osc')
clf
subplot(3,1,1)
plot(t_vec,pos_vec)
subplot(3,1,2)
plot(t_vec,deriv_vec,'k')
fun_pointer=@(t) osc_trap_pos(t,inputs);
num_deriv=num_grad_vecs(fun_pointer,t_vec,1e-5,true);
%num_deriv=derivest(fun_pointer,t_vec);
hold on
plot(t_vec,num_deriv,'r')
hold off

deriv_err=deriv_vec-num_deriv;
fprintf('rms err %G\n',rms(deriv_err))
subplot(3,1,3)
plot(t_vec,deriv_err,'k')
%%