function [F, dfM, dfE, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind, dfM_adj, dfE_adj] = teg_core_medcomp(expanded0, X0, B, Bcoder, contvec, X_to_remove)

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
    
    % Expand with within-subject coders
    tmpX = [];
    tmpX_rem = [];
    for iX = 1:size(X0, 2),
        v = X0(:, iX);
        if number_of_dims > 1,
%             v = v - mean(v);
        end;
        XM = reshape(v, size(expanded0_red));
        for iB = 1:size(Bdummy, 2),
            b = Bdummy(:, iB);
%             b = b - mean(b);
            tmp2 = (b * ones(1, size(XM, 2))) .* XM;
            tmpX = [tmpX tmp2(:)];
        end;
        for iB = 1:size(X_to_remove, 2),
            b = X_to_remove(:, iB);
%             b = b - mean(b);
            tmp2 = (b * ones(1, size(XM, 2))) .* XM;
            tmpX_rem = [tmpX_rem tmp2(:)];
        end;
    end;
    
    X = tmpX_rem;
    if ~isempty(X),
        y = expanded0_red(:) - mean(expanded0_red(:));
        b = inv(X'*X)*X'*y;
        model = X*b;
        err = y - model;
        expanded0_red = err;
    end;
    
    X0 = tmpX;
    
%     X0 = X0 - ones(size(X0, 1), 1) * mean(X0);
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
if size(transformed0, 2) > 1,
    vars = var(transformed0);
    eps_num = sum(vars) .^ 2;
    eps_den = length(vars) * sum(vars .^ 2);
    eps0 = eps_num / eps_den;
else
    eps0 = 1;
end;

nGroups = size(Bcoder, 1);
dfM = size(X0, 2);
dfM = dfM;
dfE = (nSubj - 1 - max(0, (nGroups - 1))) * number_of_dims; % Take account of matrix trickery with Bcoder
dfM_adj = eps0 * dfM;
dfE_adj = eps0 * dfE;
MSM = SSM / dfM;
MSE = SSE / dfE;
F = MSM / MSE;
p = teg_fsig(F, dfM_adj, dfE_adj);
eta = SSM / (SSM + SSE);

