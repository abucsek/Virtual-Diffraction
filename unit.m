function vec_out=unit(vec_in)

if norm(vec_in)==0
    vec_out=vec_in;
else
    vec_out=vec_in/norm(vec_in);
end