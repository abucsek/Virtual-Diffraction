function scal_out=ang(vec1_in,vec2_in)

scal_out=acos(dot(vec1_in,vec2_in)/(norm(vec1_in)*norm(vec2_in))); %(rad)
