function [pv, pv_b] = test_binom

nIts = 500000;
nVar = 4;
nVarCtrl = 4;

for N = [80], % 100 500 1000 5000 10000],
    
    pv = [];
    pv_b = [];
    for iIt = 1:nIts,
        if mod(iIt, 50) == 0,
            fprintf([num2str(N) '\t' num2str(iIt) ' / ' num2str(nIts) '\n']);
        end;
        X0 = [];
        if nVarCtrl > 0,
            X0 = randn(N, nVarCtrl);
        end;
        X0 = [ones(N, 1) X0];
        X = randn(N, nVar);
        y = floor(2 * rand(N, 1));
        [b, F, df1, df2, p, R2, pred, se_b, t_b, p_b] = teg_regression(X, y, X0);
        pv = [pv; p];
        pv_b = [pv_b; p_b(1)];
        
    end;
    
    sv = pv < 0.05;
    fprintf([num2str(N) ':\t' num2str(length(find(sv)) / length(pv)) ' (' num2str(sqrt(var(sv)) / sqrt(N)) ')\n']);
    
end;
