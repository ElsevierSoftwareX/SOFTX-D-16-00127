function [F, p, df1, df2, b_old, b_new] = teg_regr_comp(X0, Xadd, y)

f = find(isnan(mean([X0 Xadd y]'))');
X0(f, :) = [];
Xadd(f, :) = [];
y(f) = [];

[b_old, F_old, df1_old, df2_old, p_old, R2_old, pred_old] = teg_regression(X0, y, 0);

[b_new, F_new, df1_new, df2_new, p_new, R2_new, pred_new] = teg_regression([X0 Xadd], y, 0);

Err_new = var(y - pred_new);
Err_old = var(y - pred_old);
p_new = size(X0, 2) + size(Xadd, 2);
p_old = size(X0, 2);
N = length(y);
[F, p, df1, df2] = teg_ftest(Err_new, Err_old, p_new, p_old, N);
