function teg_corr_table(varargin)

% function teg_corr_table(M[, varnames])

M = varargin{1};
if length(varargin) < 2,
    varnames = {};
    for iv = 1:size(M, 2),
        varnames{iv} = ['V' num2str(iv)];
    end;
else,
    varnames = varargin{2};
end;

vc = 1:length(varnames);
nV = length(vc);
nTests = (nV * nV - nV) / 2;

corrstr = {};
for iv1 = 1:length(vc),
    for iv2 = 1:length(vc),
        M2 = M(:, [iv1 iv2]);
        f = find(isnan(mean(M2')'));
        M2(f, :) = [];
        [C, P] = corr(M2);
        if P(1, 2) < 0.05 / nTests,
            siglab = ' **';
        elseif P(1, 2) < 0.05,
            siglab = ' *';
        else,
            siglab = '';
        end;
        str0 = [num2str(C(1, 2), 2) siglab];
        corrstr{iv1, iv2} = str0;
    end;
end;

maxlen = 1;
for iv = 1:length(vc),
    l0 = length(varnames{vc(iv)});
    if l0 > maxlen,
        maxlen = l0;
    end;
end;
for iv1 = 1:length(vc),
    for iv2 = 1:length(vc),
        l0 = length(corrstr{iv1, iv2});
        if l0 > maxlen,
            maxlen = l0;
        end;
    end;
end;

eqlenlabs = {};
for iv = 1:length(vc),
    thisvn = varnames{vc(iv)};
    eqlenlabs{iv} = thisvn;
    filler0 = ' ';
    for il = 1:(maxlen - length(thisvn)),
        eqlenlabs{iv} = [eqlenlabs{iv} filler0];
    end;
end;
for iv1 = 1:length(vc),
    for iv2 = 1:length(vc),
        thiscs = corrstr{iv1, iv2};
        if mod(iv1, 2) == 0,
            filler0 = ' ';
        else,
            filler0 = ' ';
        end;
        for il = 1:(maxlen - length(thiscs)),
            corrstr{iv1, iv2} = [corrstr{iv1, iv2} filler0];
        end;
    end;
end;

for il = 1:maxlen,
    fprintf([' ']);
end;
fprintf(['\t']);
for iv1 = 1:length(vc),
    fprintf([eqlenlabs{iv1} '\t']);
end;
fprintf('\n');

for iv1 = 1:length(vc),
    for iv2 = 1:length(vc),
        if iv2 == 1,
            fprintf([eqlenlabs{iv1} '\t']);
        end;
        str0 = corrstr{iv1, iv2};
        fprintf([str0 '\t']);
    end;
    fprintf('\n');
end;
