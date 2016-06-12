% Created by Ashley Bucsek, Colorado School of Mines, 2015
% Returns a rotation matrix given a quaternion

function rotmat = quat2rot(q)
w=q(1);  x=q(2);  y=q(3);  z=q(4);
rotmat = [1-2*y^2-2*z^2  2*x*y-2*w*z  2*x*z+2*w*y;
    2*x*y+2*w*z  1-2*x^2-2*z^2  2*y*z-2*w*x;
    2*x*z-2*w*y  2*y*z+2*w*x  1-2*x^2-2*y^2];