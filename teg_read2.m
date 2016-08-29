function [D, varnames, vn2c] = teg_read2(fn)

skiplines = 0;
varnames = {};

fid = fopen(fn, 'r');
l{1} = fgetl(fid);
l{2} = fgetl(fid);
words{1} = regexp(l{1}, '\t', 'split');
words{2} = regexp(l{2}, '\t', 'split');
fclose(fid);

N = [0 0];
nNum = [0 0];
for n = 1:2,
    for m = 1:length(words{n}),
        x = str2num(words{n}{m});
        if ~isempty(x),
            nNum(n) = nNum(n) + 1;
        end;
        N(n) = N(n) + 1;
    end;
end;
r = nNum ./ N;
if r(1) < 0.5 * r(2),
    varnames = words{1};
    for ivn = 1:length(varnames),
        f = strfind(varnames{ivn}, '"');
        varnames{ivn}(f) = [];
        f = strfind(varnames{ivn}, '$');
        varnames{ivn}(f) = [];
    end;
    skiplines = 1;
end;

fid = fopen(fn, 'r');
D = [];
for n = 1:skiplines,
    l0 = fgetl(fid);
end;
while ~feof(fid),
    l0 = fgetl(fid);
    words = regexp(l0, '\t', 'split');
    tmp = [];
    for iw = 1:length(words),
        val = str2num(words{iw});
        if isempty(val),
            val = NaN;
        end;
        tmp = [tmp val];
    end;
    D = [D; tmp];
end;
fclose(fid);

vn2c = struct;
for n = 1:length(varnames),
    v = varnames{n};
    while isfield(vn2c, v),
        v = [v '_r'];
    end;
    if ~isempty(v),
        vn2c = setfield(vn2c, v, n);
    end;
end;
