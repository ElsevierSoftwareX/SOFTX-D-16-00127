function indv = teg_replace_by_order(v)

fnn = find(~isnan(v));
indv = NaN * zeros(size(v));
v = v(fnn);

[s, ind0] = sort(v);
range0 = 1:length(v);
indv(fnn(ind0)) = range0;
