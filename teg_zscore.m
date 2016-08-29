function M = teg_zscore(M)

M = M - ones(size(M, 1), 1) * mean(M);
M = ones(size(M, 1), 1) ./ sqrt(var(M));