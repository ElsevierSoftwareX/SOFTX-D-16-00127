function teg_plot_table(m, se)

% function teg_plot_table(m, se)
%
% Nesting is columns in rows.
% Rows are plotted on the main X axis.
% Columns are plotted within the clusters.
%
% That is: the values of each row will be grouped together as a cluster of bars.

if size(m, 1) == 1,
    m = m';
    se = se';
end;

xminall = Inf;
xmaxall = -Inf;
yminall = Inf;
ymaxall = -Inf;

b0 = bar(m);
hold on;
for ib = 1:length(b0),
%     c0 = get(get(b0(ib), 'Children'));
    c0 = b0(ib);
    tops0 = max(c0.YData);
    bottoms0 = min(c0.YData);
    min0x = min(c0.XData);
    if min(min0x) < xminall,
        xminall = min(min0x);
    end;
    max0x = max(c0.XData);
    if max(max0x) > xmaxall,
        xmaxall = max(max0x);
    end;
    x0v = c0.XData;
    for isb = 1:length(c0.XData),
        x0 = x0v(isb) - 0.125 + (ib - 1) * 0.5 / length(b0);
        y0 = c0.YData(isb);
        dw = 0.1;
        dh = se(isb, ib);
        plot([x0 x0], [y0 y0 + dh], 'k-');
        plot([x0 - dw x0 + dw], [y0 + dh y0 + dh], 'k-');
%         if abs(tops0(isb)) > abs(bottoms0(isb)),
%             plot([x0 x0], [y0 y0 + dh], 'k-');
%             plot([x0 - dw x0 + dw], [y0 + dh y0 + dh], 'k-');
%         else,
%             plot([x0 x0], [y0 y0 - dh], 'k-');
%             plot([x0 - dw x0 + dw], [y0 - dh y0 - dh], 'k-');
%         end;
        if y0 - dh < yminall,
            yminall = y0 - dh;
        end;
        if y0 + dh > ymaxall,
            ymaxall = y0 + dh;
        end;
    end;
end;

set(gca, 'XTick', 1:size(m, 1));

xl0 = xlim;
xlim([xl0(1) - 0.5, xl0(end) + 0.5]);
% xbuff = (xmaxall - xminall) * 0.05;
% xlim([xminall - xbuff, xmaxall + xbuff]);
ybuff = (ymaxall - yminall) * 0.05;
if ybuff > 0,
    ylim([yminall - ybuff ymaxall + ybuff]);
end;

fprintf(['print -dtiff -r600 plotfilename\n']);
