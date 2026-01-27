// The following is attached only for verifying absolutely simple
AttachSpec("../libs/sha-7-examples/src/spec");

print "**************************************************";
print "      PROOF OF THEOREM 1.4 (Magma part)           ";
print "**************************************************\n";

print "NOTE: The ideals \mathfrak{p} and \mathfrak{q} are";
print "      checked to be principal in the file ";
print "      `2-3-prop.m` because we already unpack";
print "      the ideals there. In fact, the maximal ";
print "      order in K_B has class number 1 in all cases.";

// Import the data from Table B.2
table := Read("../data/B-2-table.txt");
table := Split(table);
table := [Split(row, ":") : row in table];

for eg in table do
  // Unpack the data from one row of Table B.2
  p := eval eg[1];
  lab_A := eg[2];
  A := eval eg[3];
  dim_A := eval eg[4];
  rk_A := eval eg[5];
  cp_A := eval eg[6];
  lab_B := eg[7];
  B := eval eg[8];
  dim_B := eval eg[9];
  rk_B := eval eg[10];
  cp_B := eval eg[11];
  print "\n==================================================";
  printf "Considering the example p=%o with A=%o and B=%o\n", p, lab_A, lab_B;

  // Setup the abelian varieties (actually the curves)
  if dim_A eq 1 then
    A := EllipticCurve(A);
  else
    A := [Polynomial(a) : a in A];
    A := HyperellipticCurve(A[1], A[2]);
  end if;

  if dim_B eq 1 then
    B := EllipticCurve(B);
  else
    B := [Polynomial(b) : b in B];
    B := HyperellipticCurve(B[1], B[2]);
  end if;

  // Prove absolutely simple using criterion of Stoll (implemented in the repo
  // sha-7-examples)
  if dim_A eq 1 then
    print "- Proved A is abs. simple";
  else
    // This calls an intrinsic from `sha-7-examples`
    assert ProvedAbsolutelySimpleJacobian(A, 200);
    print "- Proved A is abs. simple";
  end if;

  // Prove the rank of A is 0
  if dim_A eq 1 then
    assert Rank(A) eq 0;
    print "- Proved the rank of A is 0";
  else
    J := Jacobian(A);
    assert RankBound(J) eq 0;
    print "- Proved the rank of A is 0";
  end if;

  // Check the torsion subgroups are trivial
  for X in <A,B> do
    if Genus(X) eq 1 then
      J := X;
    else
      J := Jacobian(SimplifiedModel(X));
    end if;
  
    assert GCD(#TorsionSubgroup(J), p) eq 1;
    print "- J[p] is trivial";
  end for;
  
  // Prove the rank of B is correct as claimed
  if Genus(B) eq 1 then
    assert Rank(B) eq rk_B;
    printf "- Proved the rank of B is %o\n", rk_B;
  else
    // The bound of 3000 is sufficient to find enough points in all cases.
    // Since J(QQ) \otimes QQ is an K-vector space it is an even rank Q-vector
    // space. So it suffices to find rk_B - 1 independent points.
    J := Jacobian(B);
    basis := ReducedBasis(Points(J : Bound:=3000));
    assert forall{P : P in basis | Order(P) eq 0}; // check infinite order
    assert #basis ge rk_B - 1;
    assert RankBound(J) eq rk_B;
    printf "- The rank of B is %o\n", rk_B;
  end if;    
end for;
