function [d, outliers] = teg_outliers(X)

for iV = 1:size(X, 2),
    vec = X(:, iV);
    m = teg_nanmean(vec);
    f = isnan(vec);
    vec(f) = m;
    X(:, iV) = vec;
end;
d = inner_get_MD(X);

ddist0 = inner_get_pdist(X);

ind = floor(0.95 * length(ddist0));
dcrit = ddist0(ind);

outliers = find(d > dcrit);

fprintf([num2str(length(outliers)) ' outliers found in N = ' num2str(size(X, 1)) '.\n']);

function ddist = inner_get_pdist(X)
nIts = 100;
n = size(X, 1);
ddist = [];
for iIt = 1:nIts,
%     selvec = 1 + floor(rand(n, 1) * n);
%     X0 = X(selvec, :);
%     d0 = inner_get_MD(X0);
%     d0crit = d0(:);
    
    X0 = randn(size(X));
    d0 = inner_get_MD(X0);
    d0 = sort(d0);
    d0crit = d0(end);
    
    ddist = [ddist; d0crit];
end;
ddist = sort(ddist);

function d = inner_get_MD(X)
[n, p] = size(X);
m0 = mean(X);
S = cov(X);

for ir = 1:n,
    d(ir) = ((X(ir, :) - m0) * inv(S) * (X(ir, :) - m0)') .^ 0.5;
end;

