function [b, F, df1, df2, p, R2, pred, se_b, t_b, p_b, R2_adj, change_R2_adj] = teg_regression(varargin)

% function [b, F, df1, df2, p, R2, pred, se_b, t_b, p_b, R2_adj, change_R2_adj] = teg_regression(X, y, X0, report, labels, b)
%
% If X0 is empty, an offset vector will be added.

X = varargin{1};
y = varargin{2};
X0 = varargin{3}; 
y_orig = y;

if length(varargin) > 3,
    report = varargin{4};
else
    report = 0;
end;

if length(varargin) > 4,
    blabels = varargin{5};
else
    blabels = {};
    for n = 1:size(X, 2),
        blabels{n} = ['b' num2str(n)];
    end;
end;

if length(varargin) > 5,
    given_b = varargin{6};
else
    given_b = [];
end;

if isempty(X0),
    X0 = ones(size(y));
end;

fisnan = find(isnan(mean([X0 X y], 2)));
fnisnan = find(~isnan(mean([X0 X y], 2)));
X(fisnan, :) = [];
y(fisnan) = [];
X0(fisnan, :) = [];
X0 = clean_X0(X0);

mean_y_orig = mean(y);
y = y - mean_y_orig;

v = var(X);
f = find(v > 0);
for iC = 1:length(f),
    X(:, f(iC)) = X(:, f(iC)) - ones(size(X, 1), 1) * mean(X(:, f(iC)));
end;

if length(unique(y)) < 2,
    b = NaN; F = NaN; df1 = NaN; df2 = NaN; p = NaN; R2 = NaN; pred = NaN; 
    se_b = NaN; t_b = NaN; p_b = NaN;
    return;
end;

if length(unique(y)) > 2,
    [b0, R20, F0, p0, t_b0, se_b0, p_b0, ErrVar0, pred0] = inner_linear(X0, y, fnisnan, mean_y_orig, y_orig, []);
    N = size(X, 1);
    k = size(X0, 2);
    if k > 1,
        k = k - 1;
        R2_adj0 = 1 - (1 - R20) * (N - 1) / (N - k - 1);
    else
        R2_adj0 = 0;
    end;
    [b, R2, F, p, t_b, p_b, se_b, ErrVar, pred] = inner_linear([X X0], y, fnisnan, mean_y_orig, y_orig, given_b);
    [F, p, df1, df2] = teg_ftest(ErrVar, ErrVar0, size(X0, 2) + size(X, 2), size(X0, 2), length(y));
    N = size(X, 1);
    k = size(X0, 2) + size(X, 2);
    if k > 1,
        k = k - 1;
        R2_adj = 1 - (1 - R2) * (N - 1) / (N - k - 1);
    else
        R2_adj = 0;
    end;
    change_R2_adj = R2_adj - R2_adj0;
    str0 = ['Change F(' num2str(df1) ', ' num2str(df2) ') = ' num2str(F) ', Change R2 = ' num2str(R2 - R20) ', Change adjusted R2 = ' num2str(change_R2_adj)];
else,
    uy = unique(y);
    biny = zeros(size(y));
    biny(y == uy(2)) = 1;
    
    if ~isempty(given_b),
        error('No given b for binomial built in.');
    end;
    [b0 dev0 stats0] = glmfit(X0, biny, 'binomial', 'link', 'logit', 'constant', 'off');
    [b, dev, stats] = glmfit([X X0], biny, 'binomial', 'link', 'logit', 'constant', 'off');
    pred = 1 ./ (1 + exp(-([X0 X] * b)));
    pred_in_orig = NaN * ones(size(y_orig));
    pred_in_orig(fnisnan) = pred;
    pred = pred_in_orig;

%     R2 = length(find(stats.resid < 0.5)) / length(biny);
    R2_adj = 0;
    change_R2_adj = 0;
    se_b = stats.se;
    F = stats.sfit;
    df1 = size(X, 2);
    df2 = (size(y, 1) - 1) - df1;
    t_b = stats.t;
    p_b = stats.p;
    
    nX = size(X, 2);
    b = b(1:nX);
    se_b = se_b(1:nX);
    t_b = t_b(1:nX);
    p_b = p_b(1:nX);
    
    Chi2stat = dev0 - dev;
    df = size(X, 2);
    p = 1 - chi2cdf(Chi2stat, df);

    R2 = Chi2stat;
    
    str0 = ['Chi2(' num2str(df) ') = ' num2str(Chi2stat) ];
end;

if report == 1,
    fprintf([str0  ', p = ' num2str(p)]); % ', R2 = ' num2str(R2)]);
    if p < 0.05,
        fprintf(' ***\n');
    else,
        fprintf(' \n');
    end;
    for ib = 1:length(blabels),
        pstr = '';
        if p_b(ib) < 0.05,
            pstr = ' *';
        end;
        fprintf(['\t' blabels{ib} '\tb = ' num2str(b(ib)) ',\tp = ' num2str(p_b(ib)) pstr '\t']);
        
        [r_solo, p_solo] = teg_corr([X(:, ib), y]);
        pstr = '';
        if p_solo(1, 2) < 0.05,
            pstr = ' *';
        end;
        fprintf(['\t(Simple corr: r = ' num2str(r_solo(1, 2)) ',\tp = ' num2str(p_solo(1, 2)) pstr ')\n']);
    end;
end;

function [b, R2, F, p, t_b, p_b, se_b, ErrVar, pred] = inner_linear(X, y, fnisnan, mean_y_orig, y_orig, given_b)
invXX = inv(X' * X);
b = invXX * X' * y;
if ~isempty(given_b),
    % mean(y) must be zero
    b(1:length(given_b)) = given_b;
    b(end) = 0;
end;
pred = X * b;
df1 = size(X, 2);
MSM = SS(pred, mean(pred)) / df1;
df2 = (size(y, 1) - 1) - df1;
MSE = SS(pred - y, 0) / df2;
F = MSM / MSE;
p = teg_fsig(F, df1, df2);
N = length(y);
se_b = sqrt(sum(y .^ 2) ./ (N - 2)) ./ sqrt(sum(X .^ 2));
se_b = se_b(:);
N = length(y);
t_b = b ./ se_b;
p_b = 2 * (1 - tcdf(abs(t_b), N - 1));
R2 = var(pred) / var(y);
ErrVar = var(y - pred);
predtmp = NaN * ones(size(y_orig));
predtmp(fnisnan) = pred;
pred = mean_y_orig + predtmp;

function X0 = clean_X0(X0)
if size(X0, 2) == 1,
    return;
end;
X0 = teg_demean(X0);
fnn = find(~isnan(mean(X0, 2)));
[O2L, Lab] = eig(cov(X0(fnn, :)));
Lab = diag(Lab);
Lab = Lab ./ mean(Lab);
fL = find(Lab >= 2 * eps);
L = X0 * O2L(:, fL);
X0 = L;

function Xc = inner_correct_X(X, X0)
Xc = [];
for iC = 1:size(X, 2),
    [b, F, df1, df2, p, R2, pred, se_b, t_b, p_b] = teg_regression(X0, X(:, iC), []);
    Xc = [Xc, X(:, iC) - pred];
end;


