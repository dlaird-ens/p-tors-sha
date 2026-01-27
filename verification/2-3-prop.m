/*
// If you are not using CHIMP (reccommended) you can uncomment this block. However
// we will not guarentee that all the dependencies are there
AttachSpec("../libs/ExactpAdics2/spec");
AttachSpec("../libs/Genus2Conductor/spec_ExactpAdics");
*/
AttachSpec("../libs/genus-2-RM/models/spec");                         //get EK model for Y_(D)
ExactpAdics_SetWarningAction("get_approx", "Ignore");                 //supress some warnings
load "modular_form_helpers.m";

////////////////////////////////////////////////////////////////////////////////
// Load the data
table_B2 := Read("../data/B-2-table_forms.txt");
table_B2 := Split(table_B2);
table_B2 := [Split(row, ":") : row in table_B2];
table_B2 := [
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
  : row in table_B2
];

////////////////////////////////////////////////////////////////////////////////
// PROOF STARTS HERE
// we first check that the modular forms really are congruent using Prop 2.2
// note that most of this code (in terms of lines) is `unpacking` the data. Most
// of the runtime is computing newspaces

for row in table_B2[#table_B2..#table_B2] do
  p := row[1]; Fp2 := GF(p^2);
  
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
      F := ResidueField(pp); Embed(F, Fp2);
      h := [Coefficient(f, ell) mod pp : ell in small_ells];
      h := [Fp2!x : x in h];
      Append(~hs_f, {[x^(p^i) : x in h] : i in [0..11]});
    end for;
  else
    pps := [p];
    h := [Fp2!Coefficient(f, ell) : ell in small_ells];
    Append(~hs_f, {h});
  end if;
  
  if Degree(Kg) ne 1 then
    qqs := [OMRepresentation(pi[1]) : pi in Factorisation(ideal<Og|p>)];
    for qq in qqs do
      F := ResidueField(qq); Embed(F, Fp2);
      h := [Coefficient(g, ell) mod qq : ell in small_ells];
      h := [Fp2!x : x in h];
      Append(~hs_g, {[x^(p^i) : x in h] : i in [0..1]});
    end for;
  else
    qqs := [p];
    h := [Fp2!Coefficient(g, ell) : ell in small_ells];
    Append(~hs_g, {h});
  end if;
  
  ////////////////////
  assert exists(ii){[i,j] : i in [1..#pps], j in [1..#qqs] |
                    #(hs_f[i] meet hs_g[j]) ge 1};
  pp := pps[ii[1]];
  qq := qqs[ii[2]];
  if Degree(Kf) ne 1 then
    assert RamificationIndex(Ideal(pp)) eq 1;
    F_pp := ResidueField(pp);
    Embed(F_pp, Fp2);
  else
    F_pp := GF(p);
  end if;
  if Degree(Kg) ne 1 then
    assert RamificationIndex(Ideal(qq)) eq 1;
    assert IsPrincipal(Ideal(qq));
    F_qq := ResidueField(qq);
    Embed(F_qq, Fp2);
  else
    F_qq := GF(p);
  end if;

  print "- The ideals pp and qq are unramified";
  print "- We also checked that the ideal qq is principal "
        cat "  (this is needed for Thm 1.4)";
  
  // Check up to sturm bound
  bound := Sturm(f, g);
  ells := PrimesInInterval(1, bound);
  flag := false;                                                      //true <=> congruence holds for all
                                                                      //    nessicary ell in the range
  for exp in [p^i : i in [0..11]] do
    flag_for_exp := true;
    for ell in ells do
      if Degree(Kf) ne 1 then
        af_modp := Fp2!(Coefficient(f, ell) mod pp);
      else
        af_modp := Fp2!Coefficient(f, ell);
      end if;
      if Degree(Kg) ne 1 then
        ag_modp := Fp2!(Coefficient(g, ell) mod qq);
      else
        ag_modp := Fp2!Coefficient(g, ell);
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
  printf "- Congruence checked for primes up to the Sturm bound %o.\n", bound;
  
  
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


////////////////////////////////////////////////////////////////////////////////
// We now check that the genus 2 curves which are claimed to be 
table_B3 := Read("../data/B-3-table.txt");
table_B3 := Split(table_B3);
table_B3 := [Split(row, ":") : row in table_B3];
table_B3 := [<row[1], row[2], eval row[3]> : row in table_B3];

EK_points := eval Read("2-3-prop_proof-data.m");

for row in table_B3 do
  ////////////////////
  N := StringToInteger(Split(row[1], ".")[1]);
  D := StringToInteger(row[2]);
  f,g := Explode([Polynomial(ff) : ff in row[3]]);
  C := HyperellipticCurve(f, g);
  i := Index(table_B3, row);
  
  print "\n==================================================";
  printf "Doing the case %o\n", row[1];
  
  // extract stored traces in other table
  assert exists(rowB2){rowB2 : rowB2 in table_B2 | row[1] in [rowB2[2], rowB2[6]]};
  if row[1] eq rowB2[2] then
    traces := rowB2[4];
  else
    traces := rowB2[8];
  end if;
  
  // indices of the good primes in PRIMES
  good_ii := [i : i in [1..30] | not PRIMES[i] in BadPrimes(C)];
  
  print "- Checking the good traces for ell < 120 agree with the stored newform";
  for i in good_ii do
    t := TraceFromCurve(C, PRIMES[i]);
    assert t eq traces[i];
  end for;
  
  print "- Checking the good traces for ell < 120 specify the newform";
  chi := DirichletCharacter(IntegerToString(N) cat ".1");
  S := CuspForms(chi, 2);
  Snew := Newforms(S);                                                //Galois orbits
  trace_agree := [];
  for f in Snew do
    if forall{i : i in good_ii | Trace(Coefficient(f[1], PRIMES[i])) eq traces[i]} then
      Append(~trace_agree, Index(Snew, f));
    end if;
  end for;
  assert #trace_agree eq 1;

  printf "- Checking the endomorphism algebras of the Jacobian\n"
         cat "  contains RM by QQ(sqrt(%o))\n", D;
  pt := EK_points[i];
  // Elkies--Kumar model is z^2 = lambda_D(g,h) check that our point satisfies this
  lambda := Getlambda_D(D);
  assert pt[3]^2 eq Evaluate(lambda, pt[1..2]);
  // Now check that the point corresponds to our curve
  I_C := IgusaClebschInvariants(C);
  I_pt := GetIgusaClebschInvariants(D : coords:=pt[1..2]);
  //IC invariants are in PP(2,4,6,10) so normalise them to affine space
  I_C_norm := [I_C[2]/I_C[1]^2, I_C[3]/I_C[1]^3, I_C[4]/I_C[1]^5];    
  I_pt_norm := [I_pt[2]/I_pt[1]^2, I_pt[3]/I_pt[1]^3, I_pt[4]/I_pt[1]^5];    
  assert I_C_norm eq I_pt_norm;
  
  ////////////////////
  print "- Checking the conductor of the Jacobian is correct";
  flag, N_J := IsSquare(Conductor_Genus2(C));
  assert flag;
  assert N_J eq N;
end for;
