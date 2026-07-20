/*
This is the same proof as in Proposition 4.4 of 
  Frengley, S. Explicit 7-torsion in the Tate–Shafarevich groups of 
  genus 2 Jacobians, J. Théor. Nombres Bordx. Volume 37 (2025) no. 2, 
  pp. 727-746. 
*/

AttachSpec("../libs/sha-7-examples/src/spec");
SetVerbose("TalkToMe", 1);

//////////////////////////////////////////////////
//               SETUP OF PROOF                 //
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// The j-invariant of the modular curve X_1(7)
function GetX1(j)
  P<t> := PolynomialRing(FieldOfFractions(Parent(j)));
  j_X1 := (t^2 + t + 1)^(3)*(t^6 - 5*t^5 - 10*t^4
           + 15*t^3 + 30*t^2 + 11*t + 1)^(3)/((t)^(7)*(t
           + 1)^(7)*(t^3 - 5*t^2 - 8*t - 1));
  return Numerator(j_X1 - j);
end function;

// The j-invariant of the modular curve X_0(7)
function GetX0(j)
  P<t> := PolynomialRing(FieldOfFractions(Parent(j)));
  j_X0 := (t^2 + 5*t + 1)^3*(t^2 + 13*t + 49)/t;
  return Numerator(j_X0 - j);
end function;

// The forgetful map X_1(7) ->  X_0(7)
function GetX1_to_X0(pt_X0)
  P<t> := PolynomialRing(FieldOfFractions(Parent(pt_X0)));
  j_X1_to_X0 := (t^3 - 5*t^2 - 8*t - 1)/(t^2 + t);
  return Numerator(j_X1_to_X0 - pt_X0);
end function;

//////////////////////////////////////////////////
//             BODY OF THE PROOF                //
//////////////////////////////////////////////////

print "**************************************************";
print "        A PROOF OF LEVEL 9510 CASE                ";
print "**************************************************\n";

// The curve C
_<xx> := PolynomialRing(Rationals());
C := HyperellipticCurve(
         Polynomial([6, -143, 1447, -5864, 11189, -10200, 3600]),
         Polynomial([-1, -1])
       );
E := EllipticCurve("9510e1");

// An alternative model for C (up to quadratic twist)
f_C := Polynomial([ 144, -528, 437, 84, 2, 12, 9 ]);
assert G2Invariants(C) eq G2Invariants(HyperellipticCurve(f_C));

// Use now the alternative model for C
C := HyperellipticCurve(f_C);
Kum := KummerSurface(Jacobian(C));
f := DefiningPolynomial(Kum);

////////////////////
// Unpack the stored proof 
proof, _ := eval Read("2-4-prop_star-cases_proof-data.m");

// unpack the fields K and L
OK := proof[1];
K := NumberField(Polynomial(OK[2]));
OK := Order([Roots(Polynomial(x), K)[1][1] : x in OK]);

// unpack the torsion of Jac(C) into min polys of the torsion point coordinates
tors_J := proof[2];
tors_J := [[OK!x : x in xi] : xi in tors_J];
tors_J := [Polynomial(g) : g in tors_J];

// the field Q(Kum[p])
L := NumberField(tors_J[3]);
Kum_L := BaseChange(Kum, L);

// unpack the torsion of Jac(C)
xi := [[r[1] : r in Roots(t, L)] : t in tors_J];
assert exists(P){Kum_L![1,xi1,xi2,xi3] : xi1 in xi[1], xi2 in xi[2], xi3 in xi[3] | Evaluate(f, [1,xi1,xi2,xi3]) eq 0};

// prove that P is 7-torsion point
assert 7*P eq Kum_L![0,0,0,1];
print "* The recorded data defines a 7-torsion point";

//////////////////// 
// Complete the proof
j := jInvariant(E);
assert ProvedSurjectiveModpGaloisRep(E, 7, 1000);
print "* The E/Q has surjective rho_7";
assert exists(r_0){r_0[1] : r_0 in Roots(GetX0(j), K)};
assert exists(r_1){r_1[1] : r_1 in Roots(GetX1_to_X0(r_0), L)};
// the following is by construction anyway
assert Evaluate(GetX1(j), r_1) eq 0;
print "* The fields Q(x(P)) are correct";
print "* This implies E[p]/{\\pm 1} \cong J[pp]/{\\pm 1}";

// Now finally check that we have the correct quadratic twist. If the quadratic
// twist is E^d (so that J[pp] \cong E^d[p]) then d is supported on bad primes 
// and 7
d := &*[a : a in BadPrimes(C) cat BadPrimes(E) cat [7]];
d := &*[a[1] : a in Factorisation(d)];
dd := Divisors(d);
dd := dd cat [-d : d in dd];

for d in dd do
  failed_ell := [];
  for ell in PrimesInInterval(13,300) do
    ap_E := TraceOfFrobenius(QuadraticTwist(E,d), ell);
    n1 := #Points(ChangeRing(C, GF(ell)));
    n2 := #Points(ChangeRing(C, GF(ell^2)));

    tr_ap_C := ell + 1 - n1;
    nm_ap_C := ExactQuotient((n1^2 + n2), 2) - (ell + 1)*n1 - ell;
    f1 := (ap_E^2 - tr_ap_C*ap_E + nm_ap_C) mod 7;
    if f1 ne 0 then
      Append(~failed_ell, ell);
    end if;
  end for;
  if d ne 1 then
    assert #failed_ell gt 0;
  end if;
end for;
print "* Checked quadratic twists supported on bad primes. \n "
      cat "    None work so this implies E[p] \cong J[pp]";



print "\n\n";
print "**************************************************";
print "        A PROOF OF LEVEL 6962 CASE                ";
print "**************************************************\n";

// The curve C
_<xx> := PolynomialRing(Rationals());
f_C := Polynomial([-708,-4484,-4071,10502,-295,-1652,-236]);
C := HyperellipticCurve(f_C);
D := HyperellipticCurve(
         Polynomial([64,176,155,34,-21,-16]), 
         Polynomial([0,-1,-1])
       );

// An alternative model for C (up to quadratic twist)
f_D := Polynomial([ 256, 704, 621, 138, -83, -64 ]);
assert G2Invariants(D) eq G2Invariants(HyperellipticCurve(f_D));

// Use now the alternative model for D
D := HyperellipticCurve(f_D);

KumC := KummerSurface(Jacobian(C));
fC := DefiningPolynomial(KumC);
KumD := KummerSurface(Jacobian(D));
fD := DefiningPolynomial(KumD);

////////////////////
// Unpack the stored proof
_, proof := eval Read("2-4-prop_star-cases_proof-data.m");

// unpack the fields K and L
OK := proof[1];
K := NumberField(Polynomial(OK[2]));
OK := Order([Roots(Polynomial(x), K)[1][1] : x in OK]);

// unpack the torsion of Jac(C) into min polys of the torsion point coordinates
tors_JC := proof[2];
tors_JC := [[OK!x : x in xi] : xi in tors_JC];
tors_JC := [Polynomial(g) : g in tors_JC];

// unpack the torsion of Jac(D) into min polys of the torsion point coordinates
tors_JD := proof[3];
tors_JD := [[OK!x : x in xi] : xi in tors_JD];
tors_JD := [Polynomial(g) : g in tors_JD];

// the field Q(KumC[p]) (which will be Q(KumD[p]))
L := NumberField(tors_JC[3]);
KumC_L := BaseChange(KumC, L);
KumD_L := BaseChange(KumD, L);

// unpack the torsion of Jac(C)
xi := [[r[1] : r in Roots(t, L)] : t in tors_JC];
assert exists(PC){KumC_L![1,xi1,xi2,xi3] : xi1 in xi[1], xi2 in xi[2], xi3 in xi[3] | Evaluate(fC, [1,xi1,xi2,xi3]) eq 0};

// unpack the torsion of Jac(D)
xi := [[r[1] : r in Roots(t, L)] : t in tors_JD];
assert exists(PD){KumD_L![1,xi1,xi2,xi3] : xi1 in xi[1], xi2 in xi[2], xi3 in xi[3] | Evaluate(fD, [1,xi1,xi2,xi3]) eq 0};


// prove that P is an 11-torsion point on KumC
assert 11*PC eq KumC_L![0,0,0,1];
print "* The recorded data defines an 11-torsion point \n"
      cat "    on Jac(C)/{\\pm 1}";
print "* We have, by construction, Q(x(P_C)) \\subset L";

// prove that P is an 11-torsion point on KumD
assert 5*PD eq 6*PD;
print "* The recorded data defines an 11-torsion point \n"
      cat "    on Jac(D)/{\\pm 1}";
print "* We have, by construction, Q(x(P_D)) \\subset L";

// It now remains to check that the Galois reps are surjective
// Let M = Jac(C)[fp] and a_ell(M) be the trace of Frob_ell 
// acting on M, for good ell. Write 
// f_ell = x^2 - a_ell(J) * x + ell.  
// Then to check surjectivity we just need to find 2 primes ell: 
// (1) such that f_ell has distinct roots mod 11
// (2) such that f_ell is irreducible mod 11
// 
// The thing to notice is that (after x -> 1/x) f_ell divides 
// the L-factor of C at ell, so we can just look for ell such 
// that L_ell(C, T) gives the hypotheses of [Ser72, Prop. 19] where
// [Ser72] J. P. Serre, Propriétés galoisiennes des points d'ordre 
//           fini des courbes elliptiques, Invent. Math. 15 (1972), 
//           no. 4, 259--331. MR387283.

function SatisfiesSerreProp19i(f)
  // f should be a quadratic polynomial over FF_p
  assert IsMonic(f);
  a := -Coefficients(f)[2];
  ell := Coefficients(f)[1];
  return IsSquare(a^2 - 4*ell);
end function;

function SatisfiesSerreProp19ii(f)
  // f should be a quadratic polynomial over FF_p
  assert IsMonic(f);
  a := -Coefficients(f)[2];
  ell := Coefficients(f)[1];
  return not IsSquare(a^2 - 4*ell);
end function;

function SatisfiesSerreProp19iii(f)
  // f should be a quadratic polynomial over FF_p
  assert IsMonic(f);
  a := -Coefficients(f)[2];
  ell := Coefficients(f)[1];
  u := a^2 / ell;
  flag := true;
  flag := flag and (not u in [0,1,2,4]);
  flag := flag and ((u^2 - 3*u + 1) ne 0);
  return flag;
end function;

function FrobPoly(X, ell)
  return EulerFactor(X, ell)^Matrix(Integers(), 2, 2, [0,1,1,0]);
end function;

P := PolynomialRing(GF(11));
ell_i := 101;
ell_ii := 271;

// 101 satisifes [Ser72, Prop 19(i)] for all possible quadratic factors
f := P!FrobPoly(C, ell_i);
ff := Factorisation(f); assert #ff eq 4; // completely split
quad_ff := [f[1]*g[1] : f,g in ff | Coefficients(f[1]*g[1])[1] eq ell_i];
assert forall{fs : fs in quad_ff | SatisfiesSerreProp19i(fs)};
// 271 satisifes [Ser72, Prop 19(ii) and (iii)] for all possible quadratic factors
f := P!FrobPoly(C, ell_ii);
ff := Factorisation(f);
assert forall{f : f in ff | Degree(f[1]) eq 2}; // product of 2 quadratics
assert forall{fs : fs in ff |
              SatisfiesSerreProp19ii(fs[1]) and
              SatisfiesSerreProp19iii(fs[1])
             };
print "* The mod pp Galois rep of Jac(C)/Q is surjective";


// 101 satisifes [Ser72, Prop 19(i)] for all possible quadratic factors
g := P!FrobPoly(D, ell_i);
gg := Factorisation(g); assert #gg eq 4; // completely split
quad_gg := [f[1]*g[1] : f,g in gg | Coefficients(f[1]*g[1])[1] eq ell_i];
assert forall{gs : gs in quad_gg | SatisfiesSerreProp19i(gs)};
// 271 satisifes [Ser72, Prop 19(ii) and (iii)] for all possible quadratic factors
g := P!FrobPoly(D, ell_ii);
gg := Factorisation(g);
assert forall{g : g in gg | Degree(g[1]) eq 2}; // product of 2 quadratics
assert forall{gs : gs in gg |
              SatisfiesSerreProp19ii(gs[1]) and
              SatisfiesSerreProp19iii(gs[1])
             };
print "* The mod pp Galois rep of Jac(D)/Q is surjective";
print "* This implies A[pp]/{\\pm 1} \cong B[qq]/{\\pm 1}";


// Now finally check that we have the correct quadratic twist. If the quadratic
// twist is C^d (so that Jac(C^d)[pp] \cong Jac(D)[qq]) 
// then d is supported on bad primes  and 11
d := &*[a : a in BadPrimes(C) cat BadPrimes(D) cat [11]];
d := &*[a[1] : a in Factorisation(d)];
dd := Divisors(d);
dd := dd cat [-d : d in dd];

for d in dd do
  failed_ell := [];
  for ell in PrimesInInterval(60,600) do
    Cd := QuadraticTwist(C, d);
    fCd := P!FrobPoly(Cd, ell);
    fCd := [fct[1] : fct in Factorisation(fCd)];
    fD := P!FrobPoly(D, ell);
    fD := [fct[1] : fct in Factorisation(fD)];
    if forall{fct : fct in fCd cat fD | Degree(fct) eq 2} then
      if not exists{fct : fct in fCd | fct in fD} then
        Append(~failed_ell, ell);
      end if;
    end if;           
  end for;
  if d ne 1 then
    assert #failed_ell gt 0;
  else
    assert #failed_ell eq 0;
  end if;
end for;
print "* Checked quadratic twists supported on bad primes. \n"
      cat "    None work so this implies A[pp] \cong B[qq]";

