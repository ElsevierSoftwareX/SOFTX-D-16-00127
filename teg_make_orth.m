function C = teg_make_orth(C)

C = demean(C);
[O2L, EV] = eig(cov(C));
C = C * O2L;
