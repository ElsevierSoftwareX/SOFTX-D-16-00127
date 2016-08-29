function teg_report(prestr, labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, p_fwe, y_red_M, fname, verbose0, bModel, X_red_M, pCritForFurther, labelFull, Btmp, Btypes, contvec)

% function teg_report(labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, p_fwe, y_red_M, fname, verbose0, bModel, X_red_M, pCritForFurther)

if verbose0 < 0,
    return;
end;

digprec1 = 3;
digprec2 = 2;
digprec3 = 3;

% resstr = [prestr fname ': <' num2str(iPred) '>' '\t' labelFull ': F(' num2str(df1, digprec) ', ' num2str(df2, digprec) ') = ' num2str(F, digprec) ', p = ' num2str(p, digprec) ', eta_p^2 = ' num2str(eta, digprec)];
resstr = [labelFull ': F(' num2str(df1, digprec1) ', ' num2str(df2, digprec1) ') = ' num2str(F, digprec2) ', p = ' num2str(p, digprec3) ',\n\teta_p^2 = ' num2str(eta, digprec2)];
if verbose0 > 0,
    resstr = [resstr '\n' prestr  '\t\tSSM = ' num2str(SSM, digprec) '\tSSE = ' num2str(SSE, digprec) '\tMSM = ' num2str(MSM, digprec) '\tMSE = ' num2str(MSE, digprec)];
    resstr = [resstr '\t epsilon = ' num2str(eps0, digprec)];
end;
fprintf([resstr]);
if p <= 0.05,
    fprintf(' * ');
end;
fprintf('\n');
fprintf([prestr '\tDescriptives:\n']);
uB = unique(Btmp);
for iGroup = 0:size(uB, 1),
    if iGroup > 0,
        fprintf(['\tGroup ' num2str(iGroup, digprec) ':\t group-code = ' num2str(Btypes(iGroup, :), digprec) '\n']);
        fg = find(Btmp == uB(iGroup));
    end;
    if iGroup == 0,
        if ~isempty(Btmp),
            continue;
        else
            fprintf('\t');
            fg = 1:size(y_red_M, 1);
        end;
    end;
    y_red_Mg = y_red_M(fg, :);
    y_red_Mw = y_red_Mg - mean(y_red_Mg, 2) * ones(1, size(y_red_Mg, 2));
    if isempty(contvec),
        f = ~isnan(mean(y_red_Mg, 2));
        myr = mean(y_red_Mg(f, :));
        if iPred > 0,
            [pvec00, t] = teg_ttest(y_red_Mw(f, :));
        else
            [pvec00, t] = teg_ttest(y_red_Mg(f, :));
        end;
%         myr = t;
        df00 = length(f) - 1;
        fprintf(['\tMeans (df = ' num2str(df00, digprec) '):\t']);
    else
        f = ~isnan(contvec(fg)) & ~isnan(mean(y_red_Mg, 2));
        [myr, pvec00] = corr(y_red_Mg(f, :), contvec(fg(f)));
        df00 = length(f) - 1;
        fprintf(['\tCorrelations (df = ' num2str(df00, digprec) '):\t']);
    end;
    for ib = 1:length(myr),
        fprintf([num2str(myr(ib), digprec)]);
        fprintf([' (p = ' num2str(pvec00(ib), digprec) ')']);
        if ib < length(myr),
            fprintf('; ');
        else
            fprintf('\n');
        end;
    end;
    if ~isempty(contvec),
        for iC = 1:size(y_red_Mg, 2),
            for iC2 = 1:(iC - 1),
                d0 = y_red_Mg(:, iC) - y_red_Mg(:, iC2);
                [myr, pvec00] = corr(d0(f), contvec(fg(f)));
                fprintf(['\t\t' num2str(iC) '-' num2str(iC2) ': r = ' num2str(myr(1, 1)) ', p = ' num2str(pvec00(1, 1), digprec) '\n']);
            end;
        end;
    else,
        for iC = 1:size(y_red_Mg, 2),
            for iC2 = 1:(iC - 1),
                d0 = y_red_Mg(:, iC) - y_red_Mg(:, iC2);
                [p, t] = teg_ttest(d0);
                fprintf(['\t\t' num2str(iC) '-' num2str(iC2) ': t(' num2str(length(d0) - 1) ') = ' num2str(t) ', p = ' num2str(p, digprec) '\n']);
            end;
        end;
    end;
end;

y_red_Mw = y_red_M - mean(y_red_M, 2) * ones(1, size(y_red_M, 2));
myr = mean(y_red_M);
seyr = var(y_red_Mw) .^ 0.5 / sqrt(size(y_red_M, 1));
if verbose0 > 0,
    fprintf([prestr '\t\tse_cells = [']);
    for ib = 1:length(myr),
        fprintf([num2str(seyr(ib), digprec)]);
        if ib < length(myr),
            fprintf(', ');
        else
            fprintf(']\n');
        end;
    end;
    tyr = myr ./ seyr;
    fprintf([prestr '\t\tt_cells = [']);
    for ib = 1:length(tyr),
        fprintf([num2str(tyr(ib), digprec)]);
        if ib < length(myr),
            fprintf(', ');
        else
            fprintf(']\n');
        end;
    end;
    fprintf([prestr '\t\tModel =[\n']);
    for ib = 1:size(X_red_M, 2),
        fprintf([prestr '\t\t']);
        fprintf([num2str(X_red_M(:, ib)')]);
        if ib < size(X_red_M, 2),
            fprintf(';\n');
        else
            fprintf(']\n');
        end;
    end;
    fprintf([prestr '\t\tb = [']);
    for ib = 1:length(bModel),
        fprintf([num2str(bModel(ib))]);
        if ib < length(bModel),
            fprintf(', ');
        else
            fprintf(']\n');
        end;
    end;
    if verbose0 > 0,
        if length(myr) > 2,
            % Paired t-tests
            fprintf([prestr '\tPost-hoc differences:\n']);
            for col1 = 1:size(y_red_M, 2),
                for col2 = 1:(col1 - 1),
                    vec = y_red_M(:, col1) - y_red_M(:, col2);
                    [p, t, df] = teg_ttest(vec);
                    if p < 0.05,
                        tresstr = [num2str(col1) ' - ' num2str(col2) ': t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)];
                        fprintf([prestr '\t\t' tresstr '\n']);
                    end;
                end;
            end;
            % vs 0
            fprintf('\tPost-hoc test vs 0:\n');
            for col1 = 1:size(y_red_M, 2),
                vec = y_red_Mw(:, col1);
                [p, t, df] = teg_ttest(vec);
                if p < 0.05,
                    tresstr = [num2str(col1) ': t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)];
                    fprintf([prestr '\t\t' tresstr '\n']);
                end;
            end;
        else,
            fprintf('\n');
        end;
        fprintf('\n');
    else,
        fprintf('\n');
    end;
else,
    fprintf('\n');
end;