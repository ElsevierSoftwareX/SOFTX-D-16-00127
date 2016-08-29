function M = demean(M)

M = M - (ones(size(M, 1), 1) * teg_nanmean(M));
