function X_to_remove = teg_X_to_remove(iBetwInt, B, bfactorStarts1, bnColsfactor1, bX1, bcellsets)

nVarsNow = length(bcellsets{iBetwInt});

factorsToRemove = [];
for n = 1:length(bfactorStarts1),
    nVars0 = length(bcellsets{n});
    if nVars0 < nVarsNow,
        factorsToRemove = [factorsToRemove n];
    end;
end;

X_to_remove = [];
for iiFactor = 1:length(factorsToRemove),
    iFactor = factorsToRemove(iiFactor);
    a = bfactorStarts1(iFactor);
    b = a + bnColsfactor1(iFactor) - 1;
    Bcoder = bX1(:, a:b);
    Bdummy = teg_B_to_BX(B, Bcoder);
    X_to_remove = [X_to_remove Bdummy];
end;
