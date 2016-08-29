function knik = teg_scree(vec)

% function knik = teg_scree(EV)
%
% EV is a vector of eigenvalues, e.g. diag(EV) where EV is from eig.

plotit = 0;

vec = vec ./ max(vec);

dvec = diff(vec);

f = find(dvec >= 1 / length(vec));
knik = f(1) + 1;
