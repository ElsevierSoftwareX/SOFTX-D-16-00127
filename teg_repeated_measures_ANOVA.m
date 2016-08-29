function O = teg_repeated_measures_ANOVA(varargin)

% function O = teg_repeated_measures_ANOVA(M, levels, varnames, Betw_group, Betw_labels, Cont, contLabels)
%
% M is observation x variable-combination matrix.
% levels is vector of levels per factor.
% varnames is a cell array of strings.
%
% Betw_group: Matrix with columns containing categorical between-subject
% factors.
% Betw_labels: cell array of strings, matching the columns of Betw-group.
%
% Betw_group: Matrix with columns containing categorical between-subject
% factors.
% Betw_labels: cell array of strings, matching the columns of Betw-group.
%
% Greenhouse-Geisser correction is applied.
%
% Thomas E. Gladwin.
% Last update: 26-08-2016.

betw_subj_interaction_depth = Inf;
pCritForFurther = 1;
pStar = 0.05;
only_show_FDR = 0;
verbose0 = 0; % -1 to remove output, 1 to provide more extensive output
plots0 = [];
digprec = 3; % digital precision of printouts.
doMANOVA = 0;

M = varargin{1};
levels = varargin{2};
pStar3 = 0.05 / prod(levels);
if size(M, 2) > prod(levels), % Number of observations per cell added
    b = size(M, 2) / 2;
    NM = M(:, (b + 1):end);
    M = M(:, 1:b);
else
    NM = ones(size(M));
end;
[nSubj, nVar] = size(M);
% Remove and store subject effects.
% Take different cell counts into account here.
subjEffect = [];
for iSubj = 1:size(M, 1),
    fNonNaN = ~isnan(M(iSubj, :));
    subjEffect(iSubj, 1) = sum(M(iSubj, fNonNaN) .* NM(iSubj, fNonNaN), 2) ./ sum(NM(iSubj, fNonNaN));
end;
M_subj = subjEffect * ones(1, size(M, 2));
M_raw = M;
M = M - M_subj;

varnames = varargin{3};
if isempty(varnames),
    varnames = cellstr(num2str((1:length(levels))'));
end;

if length(varargin) > 3,
    B = varargin{4};
    Bvarnames = varargin{5};
else,
    B = [];
    Bvarnames = {};
end;

if length(varargin) > 5,
    ContRaw = varargin{6};
    contRawLabels = varargin{7};
else,
    ContRaw = [];
    contRawLabels = {};
end;

if length(varargin) > 7,
    fname = varargin{8};
    if exist([fname '.ps']),
        delete([fname '.ps']);
    end;
else
    fname = [];
end;

if length(varargin) > 8,
    plots = varargin{9};
else
    plots = 0;
end;
if ~isempty(plots0),
    plots = plots0;
end;

if ~isempty(ContRaw),
    [Cont, contLabels, Cont_vars_involved] = teg_make_betw_cont(ContRaw, contRawLabels, betw_subj_interaction_depth);
else
    Cont = {};
end;

% Get dummy-style matrices.
[X1, factorStarts, nColsFactor, labels, cellsets, factors] = teg_create_ANOVA_dummy(levels, varnames);
y = [];
y_raw = [];
yN = [];
X = [];
for iVar = 1:nVar,
    y = [y; M(:, iVar)];
    y_raw = [y_raw; M_raw(:, iVar)];
    yN = [yN; NM(:, iVar)];
    Xsub = ones(nSubj, 1) * X1(iVar, :);
    X = [X; Xsub];
end;

O.R = [];
O.labels = {};
curr_test = [];

for iPred = 0:length(factorStarts),
    
    if iPred == 0,
        expanded0 = mean(M_raw, 2);
        y_red_M = expanded0;
        X0 = [];
        X_red_M = [];
    else,
        predcols = factorStarts(iPred):(factorStarts(iPred) + nColsFactor(iPred) - 1);
        X0 = X(:, predcols);
        X_red_M = [];
        for iX = 1:size(X0, 2),
            [X_red_M0, dum] = teg_inner_recode_raw(reshape(X0(:, iX), size(M_raw)), ones(size(M_raw)), cellsets, iPred);
            X_red_M = [X_red_M X_red_M0(1, :)'];
        end;
        
        [y_red_M, dum] = teg_inner_recode_raw(M_raw, NM, cellsets, iPred);
        M_corr = M_raw;
        subjM = teg_mean(M_corr, 2);
        M_corr = M_corr - subjM * ones(1, size(M_corr, 2));
        [dum, expanded0] = teg_inner_recode_raw(M_corr, NM, cellsets, iPred);
    end;
    
    % Between-group loop
    if isempty(B),
        bfactorStarts1 = [];
    else,
        bLevels = [];
        for iB = 1:size(B, 2);
            u = unique(B(:, iB));
            bLevels(iB) = length(u);
        end;
        [bX1, bfactorStarts1, bnColsfactor1, blabels, bcellsets, bfactors] = teg_create_ANOVA_dummy(bLevels, Bvarnames);
    end;
    for iCont = 0:size(Cont, 2),
        if iCont == 0,
            contvec = [];
            contstr = '';
        else,
            contvec = Cont(:, iCont);
            contstr = [' x ' contLabels{iCont}];
        end;
        for iBetwInt = 0:length(bfactorStarts1),
            betwStr = '';
            if iBetwInt == 0,
                Bcoder = [];
                X_to_remove = [];
            else,
                a = bfactorStarts1(iBetwInt);
                b = a + bnColsfactor1(iBetwInt) - 1;
                Bcoder = bX1(:, a:b);
                betwStr = [' x ' blabels{iBetwInt}];
                % Correct for lower-level effects and interactions
                X_to_remove = teg_X_to_remove(iBetwInt, B, bfactorStarts1, bnColsfactor1, bX1, bcellsets);
            end;
            [F, df1, df2, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind, dfM_adj, dfE_adj] = teg_RMF_ANOVA(expanded0, X0, B, Bcoder, contvec, X_to_remove);
            O.R = [O.R; F df1 df2 p MSM MSE];
            if iPred > 0,
                withinLabel = labels{iPred};
            else
                withinLabel = 'Subject-score';
            end;
            thisLabel = [num2str(size(O.R, 1)) ': ' withinLabel betwStr contstr];
            O.labels{length(O.labels) + 1} = thisLabel;
            
            if p <= pCritForFurther,
                Btmp = [];
                Btypes = [];
                if iBetwInt > 0,
                    [Btmp, Btypes] = teg_M_to_codes(B, bfactors{iBetwInt});
                end;
                teg_report('', labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, pStar3, y_red_M, fname, verbose0, bModel, X_red_M, pStar, thisLabel, Btmp, Btypes, contvec);
            end;
        end;
    end;
end;
