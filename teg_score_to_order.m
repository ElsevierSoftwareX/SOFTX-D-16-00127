function v_order = teg_score_to_order(v)

% v_order = teg_score_to_order(v)
%
% Same value -> Same order index
% NaN stays NaN

v_order = NaN * ones(size(v));

for iCol = 1:size(v, 2),
    v_order(:, iCol) = inner_func(v(:, iCol));
end;

function v_order = inner_func(v)

v_order = NaN * v;

fnnan = find(~isnan(v)); 
fnan = find(isnan(v));
v(fnan) = [];

u = unique(v);
for iu = 1:length(u),
    f = find(v == u(iu));
    v_order(fnnan(f)) = iu;
end;
