/*
This is the same proof as in Proposition 4.4 of 
  Frengley, S. Explicit 7-torsion in the Tate–Shafarevich groups of 
  genus 2 Jacobians, J. Théor. Nombres Bordx. Volume 37 (2025) no. 2, 
  pp. 727-746. 
and we verify the twists as is done in that article.
*/

AttachSpec("../libs/sha-7-examples/src/spec");
SetVerbose("TalkToMe", 1);

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
// BODY OF PROOF STARTS HERE

print "**************************************************";
print "             A PROOF OF SECTION 5                 ";
print "**************************************************\n";

// The curve C
_<xx> := PolynomialRing(Rationals());
C := HyperellipticCurve(Polynomial([-151440,91284,-26139,7086,-1113,66,-15]));

// An alternative model for C
APQ := [651/46, 23/135, 14/27];
f_C := Genus2Curve(APQ);
assert G2Invariants(C) eq G2Invariants(HyperellipticCurve(f_C));

// Use now the alternative model for C
C := HyperellipticCurve(f_C);
Kum := KummerSurface(Jacobian(C));
f := DefiningPolynomial(Kum);

// The claimed Klein quartic twists
P2<x,y,z> := ProjectiveSpace(Rationals(), 2);
Xpm := [ 
  -71*x^4 - 62*x^3*y + 42*x^3*z - 9*x^2*y^2 - 3*x^2*y*z + 42*x^2*z^2 
    + 38*x*y^3 - 39*x*y^2*z - 9*x*y*z^2 - 33*x*z^3 + 6*y^4 - 25*y^3*z + 
    6*y^2*z^2 - 13*y*z^3 + 9*z^4,
  x^4 + 11*x^3*z + 15*x^2*y^2 + 21*x^2*y*z + 9*x^2*z^2 + 21*x*y^3 
    + 21*x*y^2*z - 15*x*y*z^2 - 9*x*z^3 + 6*y^4 + 41*y^3*z 
    + 15*y^2*z^2 - 4*y*z^3 - 15*z^4
];

////////////////////
// Unpack the stored proof 
proof := eval Read("5-section_proof-data.m");

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
// Complete the proof for the twist with a point

assert DixmierOhnoInvariants(Xpm[1] : normalize:=true) eq DixmierOhnoInvariants(x^3*y + y^3*z + z^3*x : normalize:=true);
print "\n* X = X^{pm} is a twist of X(7)";
X := Curve(P2, Xpm[1]);
P := X![-1, 1/2, 1];
j := Evaluate(ModuliKleinQuarticTwist(Xpm[1]), Eltseq(P));
assert not j in [0,1728];
print "      + The point corresponds to E/Q with j \\neq 0, 1728";
assert ProvedSurjectiveModpGaloisRep(EllipticCurveWithjInvariant(j), 7, 1000);
print "   - The point corresponds to E/Q with surjective rho_7";
assert exists(r_0){r_0[1] : r_0 in Roots(GetX0(j), K)};
assert exists(r_1){r_1[1] : r_1 in Roots(GetX1_to_X0(r_0), L)};
// the following is by construction anyway
assert Evaluate(GetX1(j), r_1) eq 0;
print "   - The field Q(x(P)) is correct";

//////////////////// 
// Complete the proof for the twist without a point

assert DixmierOhnoInvariants(Xpm[2] : normalize:=true) eq DixmierOhnoInvariants(x^3*y + y^3*z + z^3*x : normalize:=true);
print "\n* X = X^{mp} is a twist of X(7)";
j_Xpm := ModuliKleinQuarticTwist(Xpm[2]);

// Build the two fields
f_1 := xx^4 - 22*xx^3 + 36*xx^2 + 72*xx - 240;
L_1 := NumberField(f_1);
f_2 := xx^4 - 22*xx^3 + 9*xx^2 - 9*xx - 224;
L_2 := NumberField(f_2);

assert not IsSquare(Discriminant(L_1)*Discriminant(L_2));
print "      + L_1 and L_2 are non-isomorphic";
assert #Subfields(L_1) eq 1; assert #Subfields(L_2) eq 1;
print "      + Fields have only Q as a subfield";

// Prove the claims for the first field
print "   - Working with field L_1";
pt := [L_1.1, 0, -2];
assert Evaluate(Xpm[2], pt) eq 0; // check it's a point

j := Evaluate(j_Xpm, pt);
assert not j in [0,1728];
print "      + The point corresponds to E/L_1 with j \\neq 0, 1728";
assert ProvedSurjectiveModpGaloisRep(EllipticCurveWithjInvariant(j), 7, 1000);
print "      + The point corresponds to E/L_1 with surjective rho_7";
  
K_L1 := ext<K | f_1>;
L_L1 := ext<L | f_1>;
Embed(L_1, K_L1, K_L1.1);
Embed(K_L1, L_L1, L_L1.1);
assert exists(r_0){r_0[1] : r_0 in Roots(GetX0(j), K_L1)};
assert exists(r_1){r_1[1] : r_1 in Roots(GetX1_to_X0(r_0), L_L1)};
// the following is by construction anyway
assert Evaluate(GetX1(j), r_1) eq 0;
print "      + The field L_1(x(P)) is correct";

// Prove the claims for the second field
print "   - Working with field L_2";
pt := [L_2.1, 1, -2];
assert Evaluate(Xpm[2], pt) eq 0; // check it's a point

j := Evaluate(j_Xpm, pt);
assert not j in [0,1728];
print "      + The point corresponds to E/L_2 with j \\neq 0, 1728";
assert ProvedSurjectiveModpGaloisRep(EllipticCurveWithjInvariant(j), 7, 1000);
print "      + The point corresponds to E/L_2 with surjective rho_7";
  
K_L2 := ext<K | f_2>;
L_L2 := ext<L | f_2>;
Embed(L_2, K_L2, K_L2.1);
Embed(K_L2, L_L2, L_L2.1);
assert exists(r_0){r_0[1] : r_0 in Roots(GetX0(j), K_L2)};
assert exists(r_1){r_1[1] : r_1 in Roots(GetX1_to_X0(r_0), L_L2)};
// the following is by construction anyway
assert Evaluate(GetX1(j), r_1) eq 0;
print "      + The field L_2(x(P)) is correct\n";
