function m = tegnanmedian(M)

if isempty(M),
    m = NaN;
    return;
end;

m = zeros(1, size(M, 2));
for iCol = 1:size(M, 2),
    vec = M(:, iCol);
    f = find(isnan(vec));
    vec(f) = [];
    m(iCol) = median(vec);
end;
