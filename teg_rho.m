function [p, rho] = teg_rho(varargin)

% [p, rho] = teg_rho(x, y)

x = varargin{1};
y = varargin{2};

[rho, p] = corr(x, y, 'type', 'Spearman');
