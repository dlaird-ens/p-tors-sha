load "modular_form_helpers.m";

////////////////////////////////////////////////////////////////////////////////
// Load the data
table := Read("../data/1-9-remark.txt");
table cat:= Read("../data/B-1-table.txt");
table := Split(table);
table := [Split(row, ":") : row in table];
table := [
  <
    eval row[1],                                                      //p
    row[2],                                                           //label of f
    eval Split(row[2], ".")[1],                                       //level of f
    eval row[3],                                                      //[ tr(a_p(f)) ]
    eval row[4],                                                      //[ nm(a_p(f)) ]
    row[5],                                                           //label of f
    eval Split(row[5], ".")[1],                                       //level of f
    eval row[6],                                                      //[ tr(a_p(f)) ]
    eval row[7]                                                       //[ nm(a_p(f)) ] 
  > 
  : row in table
];

////////////////////////////////////////////////////////////////////////////////
// Proof starts here

for row in table do
  p := row[1]; Fp12 := GF(p^12);
  
  print "\n==================================================";
  printf "Showing %o and %o are congruent (p=%o)\n", row[2], row[6], row[1];
  
  f := GetForm(row[3], row[4]);
  g := GetForm(row[7], row[8]);
  Kf := CoefficientField(f); Kf := NumberField(Kf); Of := Integers(Kf);
  Kg := CoefficientField(g); Kf := NumberField(Kf); Og := Integers(Kg);
    
  ////////////////////
  // find the congruence prime ideals
  small_ells := [ell : ell in PrimesInInterval(1, 100) | GCD(ell, row[3]*row[7]) eq 1];
  hs_f := []; hs_g := [];
  
  if Degree(Kf) ne 1 then
    pps := [OMRepresentation(pi[1]) : pi in Factorisation(ideal<Of|p>)];
    for pp in pps do
      F := ResidueField(pp); Embed(F, Fp12);
      h := [Coefficient(f, ell) mod pp : ell in small_ells];
      h := [Fp12!x : x in h];
      Append(~hs_f, {[x^(p^i) : x in h] : i in [0..11]});
    end for;
  else
    pps := [p];
    h := [Fp12!Coefficient(f, ell) : ell in small_ells];
    Append(~hs_f, {h});
  end if;
  
  if Degree(Kg) ne 1 then
    qqs := [OMRepresentation(pi[1]) : pi in Factorisation(ideal<Og|p>)];  
    for qq in qqs do
      F := ResidueField(qq); Embed(F, Fp12);
      h := [Coefficient(g, ell) mod qq : ell in small_ells];
      h := [Fp12!x : x in h];
      Append(~hs_g, {[x^(p^i) : x in h] : i in [0..11]});
    end for;  
  else
    qqs := [p];
    h := [Fp12!Coefficient(g, ell) : ell in small_ells];
    Append(~hs_g, {h});
  end if;
  
  ////////////////////
  // Proof starts here
  assert exists(ii){[i,j] : i in [1..#pps], j in [1..#qqs] | 
                    #(hs_f[i] meet hs_g[j]) ge 1};
  pp := pps[ii[1]]; 
  qq := qqs[ii[2]]; 
  if Degree(Kf) ne 1 then
    assert RamificationIndex(Ideal(pp)) eq 1;
    F_pp := ResidueField(pp); 
    Embed(F_pp, Fp12);
  else
    F_pp := GF(p);
  end if;
  if Degree(Kg) ne 1 then
    assert RamificationIndex(Ideal(qq)) eq 1;
    F_qq := ResidueField(qq);
    Embed(F_qq, Fp12);  
  else
    F_qq := GF(p);  
  end if;

  print "- The ideals pp and qq are unramified";
  
  // Check up to sturm bound
  bound := Sturm(f, g);
  ells := PrimesInInterval(1, bound);
  flag := false;                                                      //true <=> congruence holds for all 
                                                                      //    nessicary ell in the range
  for exp in [p^i : i in [0..11]] do
    flag_for_exp := true;
    for ell in ells do
      if Degree(Kf) ne 1 then
        af_modp := Fp12!(Coefficient(f, ell) mod pp);
      else
        af_modp := Fp12!Coefficient(f, ell);
      end if;
      if Degree(Kg) ne 1 then
        ag_modp := Fp12!(Coefficient(g, ell) mod qq);
      else
        ag_modp := Fp12!Coefficient(g, ell);
      end if;
      ag_modp := ag_modp^exp;                                         //apply Galois
      
      if Valuation(Level(f)*Level(g), ell) eq 0 then
        if (af_modp - ag_modp) ne 0 then
          flag_for_exp := false;
        end if;
        
      elif Valuation(Level(f)*Level(g), ell) eq 1 then
        if (af_modp*ag_modp - (ell + 1)) ne 0 then
          flag_for_exp := false;
        end if;
      end if;
    end for;
    flag := flag or flag_for_exp;                                     
  end for;
  assert flag;
  printf "- Congruence checked for primes up to the Sturm bound %o. QED\n", bound;
  
  // For fun, check Galois reps are irreducible
  ell := 2; done := false;
  _<x> := PolynomialRing(F_pp);
  while not done and (ell lt 1000) do
    if GCD(ell, Level(f)) eq 1 then
      if Degree(Kf) ne 1 then
        af_modp := F_pp!(Coefficient(f, ell) mod pp);
      else
        af_modp := GF(p)!Coefficient(f, ell);
      end if;
      done := IsIrreducible(x^2 - af_modp*x + ell);
    end if;
    ell := NextPrime(ell);
  end while;
  assert done;
  print "- Galois representation is irreducible mod pp";
end for;
