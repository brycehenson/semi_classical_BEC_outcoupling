function [pot_out,pot_grad,pot_2grad]=pot_harmonic_trap(pos,tf_param)
if size(pos,2)~=3
    error('pos must be size n x 3')
end

n_pos=size(pos,1);
square_pos_products=pos.*repmat(tf_param.inputs.omega',[n_pos,1]);
square_pos_products=square_pos_products.^2;
square_pos_products=sum(square_pos_products,2);
pot_out=(1/2)*tf_param.inputs.mass*square_pos_products;

if nargout>1
    pot_grad=pos.*(repmat(tf_param.inputs.omega',[n_pos,1]).^2);
    pot_grad=tf_param.inputs.mass*pot_grad;
end

if nargout>2
    pot_2grad=repmat(tf_param.inputs.omega',[n_pos,1]).^2;
end

end