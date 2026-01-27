/*
Magma scripts to be called from within sage. This computes the size of the 
geometric component group using RegularModel
*/
f,g := Explode([Polynomial(c) : c in C]);
C := HyperellipticCurve(f,g);
C := SimplifiedModel(C);
J := Jacobian(C);
try 
  print #ComponentGroup(RegularModel(C, p)); 
catch e 
  print -1; 
end try;
