PRIMES := PrimesInInterval(1, 120); // first 30 primes hardcoded

function eta(f, g)
  // The value eta(n) from Prop 2.2
  N := Level(f); M := Level(g);
  pp := [p[1] : p in Factorisation(N) | p[2] eq 1];
  qq := [p[1] : p in Factorisation(M) | p[2] eq 1];
  // remember have to check the sign is the same (a_ell are integers at ell|N)
  prd := &*([1] cat [p : p in qq | (p in pp) and 
                          (Integers()!Coefficient(f, p)) ne (Integers()!Coefficient(g, p))]);
  return prd * LCM(M,N);
end function;

function Sturm(f, g)
  // Makes the `Sturm` bound from Prop 2.2
  n := eta(f, g);
  mu := n*(&*[1 + 1/p[1] : p in Factorisation(n)]);
  return Floor(mu / 6);
end function;

function GetForm(level, traces)
  // Gets a weight 2 form of the right level with traces matching `traces`
  label := IntegerToString(level) cat ".1";
  chi := DirichletCharacter(label);
  S := CuspForms(chi, 2);
  Snew := Newforms(S);
  i := 1;
  found := false;
  while (not found) and (i le #Snew)  do
    ff := Snew[i];
    j := 1;
    while (not found) and (j le #ff) do
      f := ff[j];
      tt := [Trace(Coefficient(f, p)) : p in PRIMES];
      if tt eq traces then
        found := true;
        return f;
      end if;
      j +:= 1;
    end while;
    i +:= 1;
  end while;
end function;

function TraceFromCurve(C, ell)
  // Gets the value of Tr(a_l(f)) where f is the weight 2 form assocaited to Jac(C)
  F := GF(ell);
  Cell := ChangeRing(C, F);
  return ell + 1 - #Cell;
end function;
