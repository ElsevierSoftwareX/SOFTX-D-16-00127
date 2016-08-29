function [rv, pv, p] = teg_comp_corr(v1, v2, gvec, rankcor)

    % [rv, pv, p] = teg_comp_corr(v1, v2, gvec, rankcor)

    nIts = 2000;
    
    [rv, pv] = inner_c(v1, v2, gvec, rankcor);
    mad_real = inner_max_abs_diff(rv);
    
    madrv = []; 
    for iIt = 1:nIts,
        gvecRand = gvec(randperm(length(gvec)));
        [rv0, pv0] = inner_c(v1, v2, gvecRand, rankcor);
        mad_rand = inner_max_abs_diff(rv0);
        madrv = [madrv; mad_rand];
    end;

    f = find(madrv >= mad_real);
    p = length(f) / length(madrv);
    
function [rv, pv] = inner_c(v1, v2, gvec, rankcor)
    rv = [];
    pv = [];
    uG = unique(gvec);
    for iG = 1:length(uG),
        f = find(gvec == uG(iG));
        v1x = v1(f);
        v2x = v2(f);
        [C, P] = teg_corr([v1x v2x], rankcor);
        r = C(1, 2);
        p = P(1, 2);
        rv = [rv r];
        pv = [pv p];
    end;

function mad = inner_max_abs_diff(rv)
    adr = [];
    for n = 1:length(rv),
        for m = 1:(n - 1),
            adr = [adr; abs(rv(n) - rv(m))];
        end;
    end;
    mad = max(adr);
    