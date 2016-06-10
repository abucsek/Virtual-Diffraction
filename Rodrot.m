function mat_out=Rodrot(ang_in,vec_in)
% Rodrigues rotation matrix formula

v=vec_in;
theta=ang_in;

W=[0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];

mat_out=eye(3)+sin(theta)*W+2*sin(theta/2)^2*W^2;