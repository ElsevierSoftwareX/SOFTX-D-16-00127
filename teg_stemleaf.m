function O = teg_stemleaf(vecs)

% Black bar: median.
% White bar: mean.
% Red patch: range of 75% of observations around median.
% Lines: 2 * standard deviation.
% Crosses: outliers, defined as observations more than 3 * SD from the mean.

patchcol = [0.5 0.2 0.2];
if ~iscell(vecs),
    temp = vecs;
    vecs = {};
    vecs{1} = temp;
    clear temp;
end;

bw = 0.1;
bw2 = 0.75 * bw;
meanbarh = 1;
maxy = -Inf;
miny = Inf;
for n = 1:length(vecs),
    maxv = max(vecs{n});
    minv = min(vecs{n});
    if maxv > maxy,
        maxy = maxv;
    end;
    if minv < miny,
        miny = minv;
    end;
    d0 = maxv - minv;
    d0 = 0.05 * d0;
    if d0 < meanbarh,
        meanbarh = d0;
    end;
end;
mw = meanbarh;
meanbarh = mw / 2;
sdlw = bw * 0.75;
xb = 0.1;
yb = 0.1 * (maxy - miny);

all_x = [];
all_y = [];

for n = 1:length(vecs),
    xn = n;
    v0 = vecs{n};
    mv = median(v0);
    
    upv = v0(find(v0 >= mv));
    mvup = median(upv);
    upupv = upv(find(upv >= mvup));
    boxtop = median(upupv);
    
    downv = v0(find(v0 <= mv));
    mvdown = median(downv);
    downdownv = downv(find(downv <= mvdown));
    boxbottom = median(downdownv);
    
    x0 = [xn - bw, xn - bw, xn + bw, xn + bw];
    y0 = [boxbottom boxtop boxtop boxbottom];
    patch(x0, y0, patchcol);
    y0 = [mv - meanbarh mv + meanbarh mv + meanbarh mv - meanbarh];
    patch(x0, y0, [0 0 0]);
    
    O.medians(n) = mv;
    O.upper_medians(n) = boxtop;
    O.lower_medians(n) = boxbottom;
    
    mv = mean(v0);
    sd = sqrt(var(v0));
    l0 = line([xn, xn], [mv - 2 * sd, mv + 2 * sd], 'Color', [0 0 0], 'LineWidth', 2);
    l0 = line([xn - sdlw, xn + sdlw], [mv - 2 * sd, mv - 2 * sd], 'Color', [0 0 0], 'LineWidth', 2);
    l0 = line([xn - sdlw, xn + sdlw], [mv + 2 * sd, mv + 2 * sd], 'Color', [0 0 0], 'LineWidth', 2);
    all_x = [all_x; xn - sdlw; xn + sdlw];
    all_y = [all_y; mv - 2 * sd; mv + 2 * sd];
    x0 = [xn - bw2, xn - bw2, xn + bw2, xn + bw2];
    y0 = [mv - mw mv + mw mv + mw mv - mw];
    patch(x0, y0, [1 1 1]);
    
    O.means(n) = mv;
    O.sd(n) = sd;
    
    f = find(abs(v0 - mv) > 3 * sd);
    for o = 1:length(f),
        t0 = text(xn, v0(f(o)), 'X');
        set(t0, 'HorizontalAlignment', 'center');
        set(t0, 'VerticalAlignment', 'middle');
    end;
end;

set(gca, 'XTick', []);
xlabel('White = mean, black = median, box = .75 range, bars = 2 * SD');
ylim([min(all_y) - yb max(all_y) + yb]);
xlim([min(all_x) - xb max(all_x) + xb]);
