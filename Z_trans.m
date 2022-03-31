function [vec2] = Z_trans(vec)
    vec2 = log((vec+1)./(-vec+1))/2;
end