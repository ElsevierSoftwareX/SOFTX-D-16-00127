function eps0 = teg_get_eps_GG(S)

k = size(S, 1); % Number of treatments
meanVar = mean(diag(S));
matrixMean = mean(S(:));
SSMatrix = sum(S(:) .^ 2);
SSRowMeans = sum(mean(S) .^ 2);
den = k.^2 * (meanVar - matrixMean) .^ 2;
num = (k - 1) * (SSMatrix - 2 * k * SSRowMeans + k .^2 * matrixMean .^ 2);
eps0 = den / num;
