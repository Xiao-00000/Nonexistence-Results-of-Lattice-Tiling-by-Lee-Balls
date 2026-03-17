
with(combinat, partition);
with(ListTools);
with(linalg);
with(padic);
with(combinat);
with(PolynomialIdeals);
with(Groebner);
BQ := proc(r, p) local B, t, per, jiecheng, i, F, A, n, vars, BB, k, f, X, ideal, jud; B := (n, r) -> local ii; expand(add(2^ii*binomial(n, ii)*binomial(r, ii), ii = 0 .. r)); t := ceil(log[p](r + 1)); per := p^(t + 1); jiecheng := vector(p - 1, 0); jiecheng[1] := 1; for i from 2 to p - 1 do jiecheng[i] := i*jiecheng[i - 1] mod p; end do; F := X -> local i1, i2, j; expand(series((1 + q)^n*jiecheng[2*add(i1, i1 in X)]*2^numelems(X)*mul(add(i^((2*j - 1) mod (p - 1))*q^i, i = 1 .. r + 1)*x[j], j in X)/((1 - q)^(n + 1)*mul(jiecheng[2*i2], i2 in X)*mul(map2(op, 2, Statistics:-Tally(X))!~)), q, r + 1)); A := {}; for n from 0 to per - 1 do if B(n, r) mod p = 0 and B(n, r) mod p^2 <> 0 and value(coeff(F([1/2*p - 1/2]), q, r)) mod p = 0 then A := A union {n}; end if; end do; n := 'n'; if A = {} then return [r, p, per, A]; end if; vars := [seq(x[i], i = 1 .. 1/2*p - 1/2)]; BB := {}; for n in A do for k to 1/2*p - 3/2 do f[k] := add(coeff(F(X), q, r) mod p, X in partition(k)) mod p; end do; ideal := PolynomialIdeal([seq(f[i], i = 1 .. 1/2*p - 3/2)], variables = vars, characteristic = p); f[1/2*p - 1/2] := add(coeff(F(X), q, r), X in partition(1/2*p - 1/2)); try jud := timelimit(20, f[1/2*p - 1/2] in ideal); if jud then BB := BB union {n}; fprintf(fff, "%a,", n); fflush(fff); end if; catch "time expired": next; end try; end do; if numelems(BB) <> 0 then return [r, p, per, BB]; else return [r, p, per, BB]; end if; end proc;
BQ(7, 5);
fff := fopen("output.txt", WRITE);
NULL;
for i from 4 to 6 do
    for j from 10 to 20 do print(BQ(j, ithprime(i))); end do;
end do;
fclose("output.txt");

NULL;
BQ(21, 43);
evalf(23/49);
NULL;
