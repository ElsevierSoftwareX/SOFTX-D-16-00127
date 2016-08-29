function [pv, tv, dfv, effs] = teg_ttest(varargin)

% function [p, t, df, effs] = teg_ttest(vec[, report])

M = varargin{1};
report0 = 0;
if length(varargin) > 1,
    report0 = varargin{2};
end;

pv = [];
tv = [];
dfv = [];
effs = [];

for iCol = 1:size(M, 2),
    vec = M(:, iCol);
    vec(isnan(vec)) = [];
    
    t = sqrt(length(vec)) * mean(vec) / sqrt(var(vec));
    df = length(vec) - 1;
    
    p = 0.666;
    
    x = (t + sqrt(t.^2 + df)) / (2 * sqrt(t.^2 + df));
    z = df / 2;
    w = df / 2;
    try,
        tcdf00 = betainc(x, z, w);
        p = 1 - tcdf00;
    catch,
        p = NaN;
        fprintf([num2str(t) '(' num2str(df) ')\n']);
    end;
    if p > 0.5,
        p = 1 - p;
    end;
    p = 2 * p;
    
    Cohen_d = mean(vec) / sqrt(var(vec));
    
    pv = [pv p];
    tv = [tv t];
    dfv = [dfv df];
    effs = [effs Cohen_d];
    
    if report0 == 1,
        sigstr = '';
        if p < 0.05,
            sigstr = ' *';
        end
        fprintf(['t(' num2str(df) ') = ' num2str(t) ', d = ' num2str(effs) ', p = ' num2str(p) sigstr '\n']);
    end;
end;
