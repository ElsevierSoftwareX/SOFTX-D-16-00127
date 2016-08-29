function alphaSim = teg_test_ANOVA

nIts = 10000;
nSubj = 240;

alphaSim = [];
for iIt = 1:nIts,
    fprintf(['Sim ' num2str(iIt) '\n']);
    
%     levels = [3 2]; 
%     labels = {'w1', 'w3'}; 
    levels = [2 3]; 
    labels = {'A', 'B'}; 
    M = randn(nSubj, prod(levels)); 
    
    % create within-subject dependence
    for iCol = 2:size(M, 2),
        for iCol2 = 1:size(M, 2),
            if iCol2 ~= iCol,
%                 M(:, iCol) = M(:, iCol) + M(:, iCol2) / (size(M, 2) - 1);
            end;
        end;
    end;
    
%     B = [floor(4 * rand(nSubj, 1)) floor(2 * rand(nSubj, 1)) floor(3 * rand(nSubj, 1))];
%     betwLabels = {'b1', 'b2', 'b3'};
    B = [floor(4 * rand(nSubj, 2))];
    betwLabels = {'b1', 'b2'};
    
%     Cont = randn(nSubj, 3); 
%     contLabels = {'c1', 'c2', 'c3'}; 
    Cont = randn(nSubj, 2); 
    contLabels = {'c1', 'c2'}; 

    bcol = B(:, 1);
    f1 = find(bcol == 0);
    f2 = find(bcol == 1);
    wtmp1 = M(f1, 1:2:end);
    wtmp2 = M(f2, 1:2:end);
%     M(f1, 1:2:end) = M(f1, 1:2:end) + 0.5;
%     M(f2, 1:2:end) = M(f2, 1:2:end) - 0.5;
%     M(f1, 1:2:end) = wtmp1 + (Cont(f1, 2) * ones(1, size(wtmp1, 2)));
%     M(f2, 1:2:end) = wtmp2 - (Cont(f2, 2) * ones(1, size(wtmp2, 2)));

    % Create violation of sphericity
%     M(:, 2) = M(:, 1) + 4 * (M(:, 2) - M(:, 1));

    O = teg_repeated_measures_ANOVA(M, levels, labels, B, betwLabels, Cont, contLabels);
%     O = teg_repeated_measures_ANOVA(M, levels, labels, B, betwLabels);
%     O = teg_repeated_measures_ANOVA(M, levels, labels);
%     O = teg_repeated_measures_ANOVA(M, levels, labels, B, betwLabels, Cont, contLabels);
    ps = O.R(:, 4);
    alphaSim = [alphaSim; ps(:)'];
end;

disp(mean(alphaSim < 0.05, 1));
