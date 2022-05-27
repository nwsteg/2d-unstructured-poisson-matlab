% scalar multiplication of a vector
% performs vec = scalar*vec
function vec = sm(scalar,vec)
    vec.x = scalar*vec.x;
    vec.y = scalar*vec.y;
end