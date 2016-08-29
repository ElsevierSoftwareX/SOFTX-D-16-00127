function teg_plot(varargin)

% function teg_plot(C, linestyles0[, Cx, XTickLabel, title, XLabel, YLabel, legend0, remWS])
%
% Every cell of C is plotted as a line with SE-bars:
%   Let M = C{n} for some n.
%   Every row of M is an observation.
%
% Every line uses the sub-cells with parameter - value pairs in the associated cell of linestyles0.
%   E.g., linestyles0 = {{'Color', 'k', 'LineWidth', 2}, {'Color', 'r', 'LineStyle', ':', 'LineWidth', 2}};
%
% If Cx is supplied, it's either a vector of common x-values or a cell array of
% x-values per cell of C. Set Cx = [] to ignore the parameter.
%
% xticklabels: a single cell of strings corresponding to all x-values given in
% Cx, from low to high. Set as [] to ignore.
%
% title, xtitle, ytitle: strings
%
% RemWS: remove within-subject variation from calculation of error bars.
% Assumes matching matrices over all cells.

cla;

if ~iscell(varargin{1}),
    C = {varargin{1}};
else,
    C = varargin{1};
end;
linestyles0 = varargin{2};
if length(varargin) > 2,
    Cx = varargin{3};
else,
    Cx = [];
end;
if length(varargin) > 2,
    Cx = varargin{3};
else,
    Cx = [];
end;
if length(varargin) > 3,
    xticklabels = varargin{4};
else,
    xticklabels = [];
end;
if length(varargin) > 4,
    title0 = varargin{5};
else,
    title0 = [];
end;
if length(varargin) > 5,
    xlabel0 = varargin{6};
else,
    xlabel0 = 'X';
end;
if length(varargin) > 6,
    ylabel0 = varargin{7};
else,
    ylabel0 = 'Y [AU]';
end;
if length(varargin) > 7,
    legend0 = varargin{8};
else,
    legend0 = [];
end;
if length(varargin) > 8,
    RemWS = varargin{9};
else,
    RemWS = 0;
end;

inner_plots(C, linestyles0, Cx, RemWS, 0);
legend(legend0);
xvalsvec = inner_plots(C, linestyles0, Cx, RemWS, 1);

set(gca, 'XTick', unique(xvalsvec));
if ~isempty(xticklabels),
    set(gca, 'XTickLabel', xticklabels);
end;
if isempty(title0),
        title0 = ['Plot ' date ' : ' num2str(rem(now, 1))];
end;
title(title0);
xlabel(xlabel0);
ylabel(ylabel0);

function xvalsvec = inner_plots(C, linestyles0, Cx, RemWS, showBars)
C_for_se = {};
if RemWS == 1,
    subjM = zeros(size(C{1}));
    for iC = 1:length(C),
        subjM = [subjM, C{iC}];
    end;
    subjM = teg_nanmean(subjM')';
    for iC = 1:length(C),
        C_for_se{iC} = C{iC} - subjM * ones(1, size(C{iC}, 2));
    end;
else,
    C_for_se = C;
end;
SizeVals = [];
xvalsvec = [];
for iLine = 1:length(C),
    M = C{iLine};
    m = teg_nanmean(M);
    
    M_for_se = C_for_se{iLine};
    sd = sqrt(teg_nanvar(M_for_se));
    se = sd ./ sqrt(size(M_for_se, 1));
    if iscell(Cx),
        x = Cx{iLine};
    else,
        x = Cx;
    end;
    sizeVals = inner_plot(m, se, linestyles0{iLine}, x, RemWS, showBars);
    SizeVals = [SizeVals; sizeVals];
    xvalsvec = [xvalsvec; x(:)];
end;
sizeVals = [min(SizeVals(:, 1)), max(SizeVals(:, 2)), min(SizeVals(:, 3)), max(SizeVals(:, 4))];
xLims = sizeVals(1:2);
yLims = sizeVals(3:4);
dx = diff(xLims) / 10;
dy = diff(yLims) / 10;
xlim([xLims + [-dx, dx]]);
if dy > 0,
    ylim([yLims + [-dy, dy]]);
end;

function sizeVals = inner_plot(y, se, linestyles0, x, RemWS, showBars)
if isempty(x),
    x = 1:length(y);
end;
dx = mean(diff(x)) / 10;
l0 = line(x, y);
inner_set_ls(l0, linestyles0);
hold on;
if showBars == 1,
    for ix = 1:length(x),
        l0 = line([x(ix) x(ix)], [y(ix), y(ix) + se(ix)]);
        inner_set_ls(l0, linestyles0);
        l0 = line([x(ix) x(ix)], [y(ix), y(ix) - se(ix)]);
        inner_set_ls(l0, linestyles0);
        l0 = line([x(ix) - dx, x(ix) + dx], [y(ix) + se(ix), y(ix) + se(ix)]);
        inner_set_ls(l0, linestyles0);
        l0 = line([x(ix) - dx, x(ix) + dx], [y(ix) - se(ix), y(ix) - se(ix)]);
        inner_set_ls(l0, linestyles0);
    end;
end;
sizeVals = [x(1) x(end) min(y - se) max(y + se)];

function inner_set_ls(l0, linestyles0)
for ls0 = 1:2:length(linestyles0),
    set(l0, linestyles0{ls0}, linestyles0{ls0 + 1});
end;
