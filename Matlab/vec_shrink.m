function s=vec_shrink(v,a)
s=sign(v).*max((abs(v)-a),zeros(length(v), 1));
end
