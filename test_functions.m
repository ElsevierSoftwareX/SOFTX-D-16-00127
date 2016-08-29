function O = test_functions

nSubj = 150;
levels = [2 2 2];
profile0 = zeros(1, prod(levels));
for n = 1:length(levels),
    varnames{n} = ['Var' num2str(n)];
end;

profile0 = [0 1 1 0 1 0 0 1]; % effect of fac1
M = ones(nSubj, 1) * profile0;
M = M + randn(size(M));

Betw = [];
Betw_labels = {};

O = teg_repeated_measures_ANOVA(M, levels, varnames, Betw, Betw_labels);

return;

nSubj = 150;
levels = [3 3];
profile0 = zeros(1, prod(levels));
for n = 1:length(levels),
    varnames{n} = ['Var' num2str(n)];
end;
Betw = -1 + 2 * floor(rand(nSubj, 1) + 0.5);
Betw_labels{1} = 'B1';

% profile0(1:levels(end):end) = 3;
profile0 = [1 2 3 7 8 9 3 2 1]; % effect of fac1
profile0 = [0 2 0 0 8 0 0 2 0]; % effect of fac1
M = randn(nSubj, prod(levels)) + ones(nSubj, 1) * profile0;

% M = M .* (Betw * ones(1, size(M, 2)));

M = M + randn(size(M));

% Betw = [];
% Betw_labels = {};
O = teg_repeated_measures_ANOVA(M, levels, varnames, Betw, Betw_labels);
% O = teg_repeated_measures_ANOVA(M, levels, varnames);
