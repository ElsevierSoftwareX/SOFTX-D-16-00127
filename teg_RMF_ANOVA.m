function [F, dfM, dfE, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind, dfM_adj, dfE_adj] = teg_RMF_ANOVA(expanded0, X0, B, Bcoder, contvec, X_to_remove)

f = find(isnan(contvec));
if ~isempty(f),
    expanded0(f, :) = [];
    if ~isempty(X0),
        n = size(X0, 1);
        k = n / length(contvec);
        fX0 = [];
        for ik = 1:k,
            fX0 = [fX0; f + (ik - 1) * length(contvec)];
        end;
        X0(fX0, :) = [];
    end;
    contvec(f) = [];
end;

nSubj = size(expanded0, 1);
subjM = mean(expanded0, 2);

X_ind = [];
for iX = 1:size(X0, 2),
    tmp = reshape(X0(:, iX), size(expanded0));
    newcol = tmp(1, :)';
    newcol = newcol ./ sqrt(var(newcol));
    X_ind = [X_ind newcol];
end;

transformed0 = [];
number_of_dims = 1;
if ~isempty(X0),
    expanded0_red = zeros(size(expanded0));
    for iSubj = 1:size(expanded0, 1),
        y = expanded0(iSubj, :);
        y = y(:) - mean(y(:));
        
        X = X_ind;
        
        if length(find(X ~= 0)) > 0,
            b = inv(X'*X)*X'*y;
            model = X*b;
        else
            model = zeros(size(y));
        end;
        expanded0_red(iSubj, :) = model;
        
        transformed0 = [transformed0; y' * X];
        
        number_of_dims = size(X, 2);
    end;
else,
    expanded0_red = expanded0;
    X0 = ones(size(expanded0_red));
    X_ind = [];
end;

if ~isempty(Bcoder),
    Bdummy = teg_B_to_BX(B, Bcoder);
    % Adjust predictors
    nAdd = size(X0, 1) / size(Bdummy, 1);
    BdummyOrig = Bdummy;
    for iW = 2:nAdd,
        Bdummy = [Bdummy; BdummyOrig];
    end;
    Bdummy = demean(Bdummy);
    
    new_X0 = [];
    for iB = 1:size(Bdummy, 2),
        tmp = X0 .* (Bdummy(:, iB) * ones(1, size(X0, 2)));
        new_X0 = [new_X0 tmp];
    end;
    X0 = new_X0;
    
    X0 = X0 - ones(size(X0, 1), 1) * mean(X0);
end;

if ~isempty(contvec),
    contvec = contvec - mean(contvec);
    origSize = size(X0);
    tmpX = reshape(X0, nSubj, length(X0(:)) / nSubj);
    tmpX = tmpX .* (contvec * ones(1, size(tmpX, 2)));
    tmpX = reshape(tmpX, origSize);
    X0 = tmpX;
end;

X = X0;

y = expanded0_red(:) - mean(expanded0_red(:));
if isnan(rcond(X' * X)) || min(var(X)) == 0,
    b = zeros(size(X, 2), 1);
else,
    b = inv(X'*X)*X'*y;
end;
model = X*b;
err = y - model;
SSM = sum(model.^2);
SSE = sum(err.^2);

bModel = b;

% Calc eps.
if size(expanded0_red, 2) > 1,
    C = cov(expanded0_red);
    [O2L, L] = eig(C);
    L = diag(L);
    d = length(find(L > 2 * eps));
    eps0 = (sum(L) .^ 2) / (d * sum(L .^ 2));
else
    eps0 = 1;
end;

nGroups = size(Bcoder, 2) + 1;
dfM = size(X0, 2);
dfM = dfM;
dfE = (nSubj - 1 - max(0, (nGroups - 1))) * number_of_dims; % Take account of matrix trickery with Bcoder
dfM_adj = eps0 * dfM;
dfE_adj = eps0 * dfE;
MSM = SSM / dfM_adj;
MSE = SSE / dfE_adj;
F = MSM / MSE;
p = teg_fsig(F, dfM_adj, dfE_adj);
eta = SSM / (SSM + SSE);

