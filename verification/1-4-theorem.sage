from sage_cluster_pictures.cluster_pictures import Cluster
load("cp_helpers.sage")  
# Will use cluster pictures when l || N is odd 
# Will call Magma otherwise (this is what cp_helpers.sage does)

print("**************************************************")
print("      PROOF OF THEOREM 1.4 (Sage part)           ")
print("**************************************************\n")

if __name__ == "__main__":
    print("We verify the Tamagawa numbers computed in Table B.2\n")

    # Unpack the table into strings
    with open("../data/B-2-table.txt", "r") as f:
        table = f.readlines()
        table = [row.split(":") for row in table]

    for eg in table:
        # extract from the string the data in the table
        p = int(eg[0])
        lab_A = eg[1]
        A = eval(eg[2])
        dim_A = int(eg[3])
        rk_A = int(eg[4])
        cp_A = int(eg[5])                                   #= -1 if conjectural
        lab_B = eg[6]
        B = eval(eg[7])
        dim_B = int(eg[8])
        rk_B = int(eg[9])
        cp_B = int(eg[10])                                  #= -1 if conjectural

        for cp_info in [(lab_A, A, dim_A, cp_A), (lab_B, B, dim_B, cp_B)]:
            lab_J, J, dim_J, cp_J = cp_info
            print("==================================================")
            print(f"Example with label {lab_J}")
            if dim_J == 1:
                c_prod = EllipticCurve(J).tamagawa_product()
                assert cp_J == c_prod
                print(f"- The Tamagawa product is {c_prod}, as required\n")

            elif cp_J != -1:
                # Dimension 2 and we claim to compute prod(cp) (only at
                # multiplicitive primes if p >= 7)
                N = int(lab_J.split(".")[0])                # Level

                # For multiplicative primes not 2 use the cluster pictures
                mult_ps = mult_primes(N)
                c_prod = prod([clusters_component_gp(J, ell)
                               for ell in mult_ps if ell != 2])

                # If p = 2 is multiplicative, use regular model
                if 2 in mult_ps:
                    c_prod *= magma_get_component_gp_regmod(J, 2)

                # If p = 5 we must do additive primes too
                if p == 5:              
                    bad_ps = pdiv(N)
                    add_ps = [ell for ell in bad_ps if not ell in mult_ps]
                    c_prod *= prod([magma_get_component_gp_regmod(J, ell)
                                    for ell in add_ps])

                # Conclude
                assert c_prod == cp_J
                print(f"- The Tamagawa product (except at additive primes\n" +
                      f"  if p > 5) is equal to {c_prod}, as required\n")

            else:
                # Dimension 2 but we make no claim of proof
                print("- In this case we do not claim to compute the\n" + 
                      "  Tamagawa product\n")
