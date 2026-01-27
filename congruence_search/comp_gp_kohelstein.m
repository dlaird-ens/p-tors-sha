chi := IntegerToString(level) cat ".1";
chi := DirichletCharacter(chi);
S := CuspForms(chi, 2);
N := Newforms(S);
N := [fs : fs in N | Degree(BaseRing(fs[1])) eq dim];
K := NumberField(Polynomial(field_poly));
if dim ne 1 then
  N := [fs : fs in N |
        IsIsomorphic(NumberField(DefiningPolynomial(BaseRing(fs[1]))), K)
       ];
end if;
small_pp := PrimesInInterval(1, 10000)[1..#my_traces];

for fs in N do
  for f in fs do
    tr := [Trace(a) : a in [Coefficient(N[1][1], p) : p in small_pp]];
    nm := [Norm(a) : a in [Coefficient(N[1][1], p) : p in small_pp]];
    if tr eq my_traces then
      if nm eq my_norms then
        my_f := f;
      end if;
    end if;
  end for;
end for;

A := ModularAbelianVariety(my_f);
ret := [];
for p in ps do
  try
    Append(~ret, ComponentGroupOrder(A, p));
  catch e
    Append(~ret, -1);
  end try;
end for;

print ret;
