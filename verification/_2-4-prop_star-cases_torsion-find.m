/***************************************************
A script to compute pp | 11 division polynomials via
a number of tricks. 
***************************************************/

// You need to put this attach command in, see the README.md
// Attach("/path/to/tricks.m");
print "You need to attach the tricks file..., see the README.md";

function Times3Kummer(K)
  // A function to compute polynomials for the
  // multiplication by 3 map on a Kummer surface
  P3 := ProjectiveSpace(BaseRing(K), 3);

  two_times := [Evaluate(p, [P3.i : i in [1..4]]) : p in K`Delta];

  BB := K`BBMatrix;
  L12 := two_times cat [P3.i : i in [1..4]];
  L3 := [P3.i : i in [1..4]];
  c1 := P3.1; c2 := Evaluate(BB[1,1], L12);
  three_times := [ c1*Evaluate(BB[1,j], L12) - L3[j]*c2 : j in [1..4] ];
  three_times[1] := c1*c2;

  return three_times;
end function;

//////////////////////////////////////////////////
// Setup 
R<t> := PolynomialRing(Rationals());
F := Polynomial([ 256, 704, 621, 138, -83, -64 ]);
C := HyperellipticCurve(F);
Kum := KummerSurface(Jacobian(C));

CC := ComplexField(500);                      // In many cases this is not enough prec
P<x> := PolynomialRing(CC);
F_C := P!F;

// Some auxiliary data for the Kummer surface
K_ply := DefiningPolynomial(Kum);
F0 := Kum`AuxPolynomials[1];

//////////////////////////////////////////////////
// Analytic computation of the endomorphism ring

A := AnalyticJacobian(F_C);
E := EndomorphismRing(A); assert IsCommutative(E);
gg := [g : g in Generators(E)];
M := gg[1]; /* M := E!1 - gg[1];  // Might need to use Galois conjugate*/
assert Eltseq(MinimalPolynomial(M)) eq [-1, -1, 1];        // Golden ratio
period_mat := BigPeriodMatrix(A);

// Compute the action of phi on the period matrix
alpha := Submatrix(period_mat*Matrix(CC,M),1,1,2,2)*Submatrix(period_mat,1,1,2,2)^-1 ;
alpha := Matrix(CC,2,2,[Round(Real(x)) : x in Eltseq(alpha)]);


//////////////////////////////////////////////////
// Computing the images of points under the end
// of the Analytic Jacobian

// Generate some random x-coordinates on our curve C
xxQ := [[Random(-200,200)/Random(1,200) : _ in [1..2]] : _ in [1..1000]];
// Now lift them to some points on C
xx := [[
          <CC!myx[1], Roots(x^2 - Evaluate(F_C, myx[1]))[1][1]>,
          <CC!myx[2], Roots(x^2 - Evaluate(F_C, myx[2]))[1][1]>
        ]
       : myx in xxQ];
// D = P_1 + P_2 - \pi^*(infty) is an element of Jac(C) and so we write down its
// image on the Kummer surface. This gives our source points.
src := [
  [
    1,
    x[1][1] + x[2][1],
    x[1][1]*x[2][1],
    (Evaluate(F0, [x[1][1] + x[2][1], x[1][1]*x[2][1]]) - 2*x[1][2]*x[2][2])/(x[1][1] - x[2][1])^2
  ]
  : x in xx
];

// Compute the image of our pairs of points `xx` in the AnalyticJacobian `A`
// Note we get pairs for P1 - (weierstrass pt) and P2 - (weierstrass pt) so we
// need to sum them to get P1 + P2 - \pi^*(\infty)
pp := [
  <
    ToAnalyticJacobian(p[1][1], p[1][2], A),
    ToAnalyticJacobian(p[2][1], p[2][2], A)
  >
  : p in xx
];
target := [FromAnalyticJacobian(alpha*(p[1] + p[2]), A) : p in pp];
// Compute the images of these points on the Kummer surface.
target := [
  [
    1,
    s[1][1] + s[2][1],
    s[1][1]*s[2][1],
    (Evaluate(F0, [s[1][1] + s[2][1], s[1][1]*s[2][1]]) - 2*s[1][2]*s[2][2])/(s[1][1] - s[2][1])^2
  ]
  : s in target
];


//////////////////////////////////////////////////
// Setting up the points over a finite field

/*
DISCUSSION: Naively one would like to interpolate now this map just by working
over CC and approximating the coefficients by rational numbers. Unfortunately
since the resulting linear algebra is so dense, this completely fails. Instead
we have to work a little harder and work out what these points are over a number
field. Then we reduce them modulo p. The important thing is to make sure we
actually know the signs (otherwise the interpolation won't work). So even though
all these points are defined over a quadratic field, we choose p and then only
keep those points which become F_p rational.
*/


discs := [];                              // This is discs of the field defining the point
rat_src := [];                            // These will be the source points
rat_tar := [];                            // These will be the target points

// This loop will now store the discriminants of the quadratic fields and
// points as [1, a1 + a2*t, b1 + b2*t, c1 + c2*t] where t^2 = D
for i in [1..#src] do
  D := Discriminant(
           QuadraticField(Discriminant(MinimalPolynomial(src[i][4], 2)))
         );
  K := NumberField(t^2 - D);
  Append(~discs, D);
  assert K.1^2 eq D;
  
  rat_s := [];
  for s in src[i] do
    // choose the correct conjugate
    assert exists(r){
                   r[1] : r in Roots(MinimalPolynomial(s, 2), K) |
                   Abs(CC!Conjugates(r[1])[1] - s) lt 10^(-10)
                 };
    Append(~rat_s, Eltseq(r));
  end for;
  rat_t := [];
  for t in target[i] do
    // choose the correct conjugate
    assert exists(r){
                   r[1] : r in Roots(MinimalPolynomial(t, 2), K) |
                   Abs(CC!Conjugates(r[1])[1] - t) lt 10^(-10)
                 };
    Append(~rat_t, Eltseq(r));
  end for;
  Append(~rat_src, rat_s);
  Append(~rat_tar, rat_t);
end for;

// Now we coerce these points to live over F_p
prime := NextPrime(10^50);                                  // moderate prime
Fp := GF(prime);
_<tt> := PolynomialRing(Fp);
ii := [i : i in [1..#discs] | IsSquare(Fp!discs[i])];
discs := discs[ii];
rat_src := rat_src[ii];
rat_tar := rat_tar[ii];

fp_src := [];
fp_tar := [];
for i in [1..#discs] do
  d := Roots(tt^2 - discs[i])[1][1];
  fp_s := [];
  for s in rat_src[i] do
    Append(~fp_s, s[1] + d*(Fp!s[2]));
  end for;
  fp_t := [];
  for t in rat_tar[i] do
    Append(~fp_t, t[1] + d*(Fp!t[2]));
  end for;
  Append(~fp_src, fp_s);
  Append(~fp_tar, fp_t);
end for;

//////////////////////////////////////////////////
// Now the interpolation step
P<[x]> := PolynomialRing(Rationals(), 4);
F := FieldOfFractions(P);

// phi will be the equation for the endomorphism on the Kummer surface
phi := [];
for i in [2,3,4] do
  // Interpolate the map K->K->A1 where the last map is projection onto
  // the i-th coordinate. It turns out the degree is 3 but I couldn't have
  // guessed that, I just tried 2 and then 3. If you get nonsense you probably
  // need to increase `prime` or use more points.
  f := InterpolateMapToA1(fp_src, [t[i] : t in fp_tar], 3)[1];
  // Now lift the rational function f[1]/f[2] up to characteristic 0 using the
  // LLL algorithm
  f1,d := LiftPolyViaLLL(f[1]);
  f2,d2 := LiftPolyViaLLL(d*f[2]);
  f := d2*(F!f1/F!f2);
  // Now store it
  Append(~phi, f);
end for;

/////////////////////////////////////////////////
// Now the calculation of the division polynomials

// The point is that Nm(3 - phi) = 11, so we just need to compute
// an equation for [3]P = phi(P)
three := Times3Kummer(Kum);
three := [F | three[2]/three[1], three[3]/three[1], three[4]/three[1]];

fs := [P!Numerator(three[i] - phi[i]) : i in [1..3]];
fK := P!DefiningPolynomial(Kum);

A3<x,y,z> := PolynomialRing(Rationals(), 3);  _<X> := PolynomialRing(Rationals());
fs := [Evaluate(ff, [1,x,y,z]) : ff in fs];
fK := Evaluate(fK, [1,x,y,z]);

r1 := Resultant(fK, fs[1], 3);
r2 := Resultant(fK, fs[2], 3);

assert exists(R){R[1] : R in Factorisation(Evaluate(Resultant(r1, r2, 2), [X,0,0])) | Degree(R[1]) eq 60};
assert exists(S){S[1] : S in Factorisation(Evaluate(Resultant(r1, r2, 1), [0,X,0])) | Degree(S[1]) eq 60};

R := PolynomialRing(Integers())! (R*LCM([Denominator(c) : c in Coefficients(R)]));
S := PolynomialRing(Integers())! (S*LCM([Denominator(c) : c in Coefficients(S)]));


//////////////////////////////////////////////////
// Making a better presentation of this mess

// Now we have minimal polynomials for the first 2 coordinates, but unfortunately
// we still have work to do, because Magma will refuse to factor these for us
// over the fields they define (or it just takes forever). The point is now to
// compute the degree 12 = 11 + 1 subfield (the "rational cyclic subgroup" field
// like X_0(11) vs X_1(11)).

// The idea is to go back to the analytic Jacobian
M := 3*(E!1) + M;
assert Determinant(M) eq 121;
norm_11_elt := Submatrix(period_mat*Matrix(CC,M),1,1,2,2)*Submatrix(period_mat,1,1,2,2)^-1 ;
norm_11_elt := Matrix(CC,2,2,[Round(Real(x)) : x in Eltseq(norm_11_elt)]);

// Now compute a point in the kernel of the norm_11_elt
p := ColumnSubmatrixRange(period_mat, 1, 1);
p := norm_11_elt^(-1)*p;

// Now compute its span under (Z/11 Z)* / {\pm 1}
all_p := [i*p : i in [1..5]];

// Now compute their image on the kummer surface
all_pp := [FromAnalyticJacobian(p,A) : p in all_p];
xpx := [pp[1][1] + pp[2][1] : pp in all_pp];
xx := [pp[1][1] * pp[2][1] : pp in all_pp];

/*
// Now, at this points you could try to compute the R and S, and you might
// wonder why we did it our way. Unfortunately we just don't have enough complex
// precision to survive this calculation.
ply := MinimalPolynomial(xpx[1], 60);
ply2 := MinimalPolynomial(xx[1], 60);
*/

// We do have enough precision however to compute the minimal polynomials of the
// average over the orbit of the (Z/ 11)* / {\pm 1}
ply := MinimalPolynomial(&+xpx, 12);
ply2 := MinimalPolynomial(&+xx, 12);

// Make monic polynomials defining the same field
ply := ply^Matrix(Integers(), 2, 2, [1,0,0,LeadingCoefficient(ply)]);
ply2 := ply2^Matrix(Integers(), 2, 2, [1,0,0,LeadingCoefficient(ply2)]);

// Becuase we have 2 (more or less unrelated) polynomials we can compute a
// good order much faster by using them both, this quickly cleans up the
// equation for the field (otherwise OtimisedRepresentation would take ages)
F := NumberField(ply);
r1 := Roots(ply, F)[1][1];
r2 := Roots(ply2, F)[1][1];
O := Order([r1^i : i in [0..11]] cat
           [r2^i : i in [0..11]]);
O := LLL(O);
Omax := MaximalOrder(O);
Omax := LLL(Omax);

B := Basis(Omax);
B := [MinimalPolynomial(b) : b in B];
B := [Eltseq(b) : b in B];
cB := B;

//////////////////////////////////////////////////
// Now back to computing the torsion point

// Unpack B again 
f := Polynomial(B[2]);
KK := NumberField(f);
B := [Roots(Polynomial(d), KK)[1][1] : d in B];
O := Order(B); assert IsMaximal(O);
P := PolynomialRing(KK);

// Factorise R over O
R := Factorisation(P!R)[1][1];
d := &cat[Eltseq(e) : e in Coefficients(R)];
d := [Denominator(e) : e in d];
d := LCM(d);
R := R*d;
cR := [O | c : c in Coefficients(R)];

// Factorise S over O
S := Factorisation(P!S)[1][1];
d := &cat[Eltseq(e) : e in Coefficients(S)];
d := [Denominator(e) : e in d];
d := LCM(d);
S := S*d;
cS := [O | c : c in Coefficients(S)];

// Compute the final coordinate over L = QQ(x(P))
L := NumberField(R); _<tt> := PolynomialRing(L);
rr := [r[1] : r in Roots(R, L)];
ss := [s[1] : s in Roots(S, L)];
KumL := BaseChange(Kum, L);

// One of the (r,s) pairs is valid and the others are nonsense
assert exists(i){i : i in [1..5] | 
                 IsSquare(Discriminant(Evaluate(fK, [rr[1], ss[i], tt])))
                };

f := Evaluate(fK, [rr[1], ss[i], tt]);
_, sqrtd := IsSquare(Discriminant(f));
z := (-Eltseq(f)[2] + sqrtd) / (2*Eltseq(f)[3]); assert Evaluate(f, z) eq 0;
if 11*KumL![1, rr[1], ss[i], z] ne KumL![0,0,0,1] then
  z := (-Eltseq(f)[2] - sqrtd) / (2*Eltseq(f)[3]); assert Evaluate(f, z) eq 0;
end if;

T := MinimalPolynomial(z);
d := &cat[Eltseq(e) : e in Coefficients(T)];
d := [Denominator(e) : e in d];
d := LCM(d);
T := T*d;
cT := [O | c : c in Coefficients(T)];


// Now the data you want is 
print cB;
print [cR, cS, cT];
