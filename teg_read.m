function D = teg_read(fn, skiplines)

fid = fopen(fn, 'r');
D = [];
for n = 1:skiplines,
    l0 = fgetl(fid);
end;
while ~feof(fid),
    l0 = fgetl(fid);
    fs = find(isspace(l0));
    tmp = [];
    a = 1;
    for is = 1:(length(fs) + 1),
        if is <= length(fs),
            b = fs(is);
        else,
            b = length(l0) + 1;
        end;
        if a < b,
            num0 = str2num(l0(a:(b - 1)));
        else,
            num0 = NaN;
        end;
        if isempty(num0),
            num0 = NaN;
        end;
        a = b + 1;
        tmp = [tmp num0];
    end;
    e0 = length(tmp) - size(D, 2);
    if e0 > 0,
        F0 = NaN * zeros(size(D, 1), e0);
        D = [D F0];
    end;
    if e0 < 0,
        tmp = [tmp NaN * zeros(1, abs(e0))];
    end;
    D = [D; tmp];
end;
fclose(fid);
