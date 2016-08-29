function [Btmp, Btypes] = teg_M_to_codes(B, bfactors)

Btmp = B(:, bfactors);

v = [];
for n = 1:size(Btmp, 2),
    v0 = 10 .^ (n - 1);
    v = [v; v0];
end;

Btmp = Btmp * v;

u = unique(Btmp);
Btypes = [];
for iu = 1:length(u),
    f = find(Btmp == u(iu));
    Btypes = [Btypes; Btmp(f(1), :)];
end;
