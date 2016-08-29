function Bdummy = teg_B_to_BX(B, Bcoder)

% recode B to index values
Blevels = [];
for iB = 1:size(B, 2),
    bv = B(:, iB);
    u = unique(bv);
    tmp = zeros(size(bv));
    for iu = 1:length(u),
        f = find(bv == u(iu));
        tmp(f) = iu;
    end;
    B(:, iB) = tmp;
    Blevels = [Blevels length(u)];
end;
Bdummy = [];
for iSubj = 1:size(B, 1),
    grVars = B(iSubj, :);
    % Find dummy coder for this combination of group-variables
    row = 1;
    for iB = 1:size(B, 2),
        if iB < length(Blevels),
            skipper = prod(Blevels((iB + 1):end));
        else,
            skipper = 1;
        end;
        row = row + (grVars(iB) - 1) * skipper;
    end;
    newrow = Bcoder(row, :);
    Bdummy = [Bdummy; newrow];
end;
