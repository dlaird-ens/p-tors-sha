function ReturnCongs(filename)
  code := "import json\n"
          cat Sprintf("with open('%o', 'r') as f:\n", filename)
          cat "    x = json.load(f) \n"
          cat "x = [[a['label'] for a in xx] for xx in x]; print(x)";
  my_str := Pipe(Sprintf("python3 -c \"%o\"", code), "");
  
  ret := "";
  for i in [1..#my_str] do
    if my_str[i] eq "'" then
      ret cat:= "\"";
    else
      ret cat:= my_str[i];
    end if;
  end for;
  return eval ret;
end function;

function LoadWeier()
  my_str := Read("data/label_to_curve.txt");
  my_str := Split(my_str, "\n");
  my_str := [Split(s, ":") : s in my_str];
  my_str := [<s[1], eval s[2]> : s in my_str];
  lab := [s[1] : s in my_str];
  weier := [s[2] : s in my_str];
  return lab, weier;
end function;

function LabToCurve(lab, labs, weiers)
  if lab in labs then
    i := Index(labs, lab);
    return weiers[i];
  end if;
  return [];
end function;

SetClassGroupBounds("GRH");
p := 31;

labs, weiers := LoadWeier();
file := Sprintf("congruences/0/%o.json", p);
congs := ReturnCongs(file);
congs := [[<l, LabToCurve(l, labs, weiers)> : l in c] : c in congs];
congs := [c : c in congs | #[a : a in c | #a[2] ge 1] gt 1];

small_fct := [-1,2,-2,5,-5,7,-7,19,-19];

for c in congs do
  C1 := HyperellipticCurve(Polynomial(c[1][2][1]), Polynomial(c[1][2][2]));
  C2 := HyperellipticCurve(Polynomial(c[2][2][1]), Polynomial(c[2][2][2]));
  
  for d in small_fct do
    J1 := Jacobian(QuadraticTwist(C1, d));
    J2 := Jacobian(QuadraticTwist(C2, d));
    _,_,_,_,r1 := MordellWeilGroupGenus2(J1 : RankOnly:=true);
    _,_,_,_,r2 := MordellWeilGroupGenus2(J2 : RankOnly:=true);
    if r1 gt r2 then
      c[1][1], c[2][1], d;
      MordellWeilGroupGenus2(J1 : RankOnly:=true);
      MordellWeilGroupGenus2(J2 : RankOnly:=true);
    end if;
  end for;
end for;
