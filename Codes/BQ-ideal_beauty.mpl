
with(combinat, partition);

with(ListTools);

with(linalg);
with(combinat);
with(PolynomialIdeals);
with(Groebner);
NULL;
lattice_tiling_check := proc(r, p) local compute_lee_ball, per, factorial_mod, index, F, valid_ns, temp_n, vars, i, confirmed_ns, k, f, ideal, ii, X; compute_lee_ball := (n, r) -> expand(add(2^ii*binomial(n, ii)*binomial(r, ii), ii = 0 .. r)); per := p^(1 + ceil(log[p](r + 1))); factorial_mod := vector(p - 1, 0); factorial_mod[1] := 1; for index from 2 to p - 1 do factorial_mod[index] := index*factorial_mod[index - 1] mod p; end do; F := X -> local i1, i2, j; expand(series((1 + q)^temp_n*factorial_mod[2*add(i1, i1 in X)]*2^numelems(X)*mul(add(i^((2*j - 1) mod (p - 1))*q^i, i = 1 .. r + 1)*x[j], j in X)/((1 - q)^(temp_n + 1)*mul(factorial_mod[2*i2], i2 in X)*mul(map2(op, 2, Statistics:-Tally(X))!~)), q, r + 1)); valid_ns := {}; for temp_n from 0 to per - 1 do if compute_lee_ball(temp_n, r) mod p = 0 and compute_lee_ball(temp_n, r) mod p^2 <> 0 then valid_ns := valid_ns union {temp_n}; end if; end do; temp_n := 'temp_n'; vars := [seq(x[i], i = 1 .. 1/2*p - 1/2)]; confirmed_ns := {}; for temp_n in valid_ns do for k to 1/2*p - 3/2 do f[k] := add(coeff(F(X), q, r), X in partition(k)); end do; ideal := PolynomialIdeal([seq(f[i], i = 1 .. 1/2*p - 3/2)], variables = vars, characteristic = p); f[1/2*p - 1/2] := add(coeff(F(X), q, r), X in partition(1/2*p - 1/2)); if f[1/2*p - 1/2] in ideal then confirmed_ns := confirmed_ns union {temp_n}; end if; end do; return [r, p, per, confirmed_ns]; end proc;
# 
NULL;
lattice_tiling_check(7, 5);
NULL;
