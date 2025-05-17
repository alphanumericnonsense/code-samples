library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
-- test vector with intermediate values
-- https://csrc.nist.gov/CSRC/media/Projects/Cryptographic-Standards-and-Guidelines/documents/examples/SHA256.pdf
package sha2_256_pkg is
    type slv_1_t is array(natural range <>) of std_logic_vector;
    -- initial hash constants
    constant SHA2_256_INIT_ARRAY_c : slv_1_t(0 to 8-1)(0 to 32-1) :=
        (X"6a09e667", X"bb67ae85", X"3c6ef372", X"a54ff53a",
        X"510e527f", X"9b05688c", X"1f83d9ab", X"5be0cd19");
    constant SHA2_256_INIT_c : std_logic_vector(0 to 256-1) := (X"6a09e667bb67ae853c6ef372a54ff53a510e527f9b05688c1f83d9ab5be0cd19");
    -- 64 round constants
    constant SHA2_256_ROUND_c : slv_1_t(0 to 64-1)(0 to 32-1) :=
        (X"428a2f98", X"71374491", X"b5c0fbcf", X"e9b5dba5", X"3956c25b", X"59f111f1", X"923f82a4", X"ab1c5ed5",
        X"d807aa98", X"12835b01", X"243185be", X"550c7dc3", X"72be5d74", X"80deb1fe", X"9bdc06a7", X"c19bf174",
        X"e49b69c1", X"efbe4786", X"0fc19dc6", X"240ca1cc", X"2de92c6f", X"4a7484aa", X"5cb0a9dc", X"76f988da",
        X"983e5152", X"a831c66d", X"b00327c8", X"bf597fc7", X"c6e00bf3", X"d5a79147", X"06ca6351", X"14292967",
        X"27b70a85", X"2e1b2138", X"4d2c6dfc", X"53380d13", X"650a7354", X"766a0abb", X"81c2c92e", X"92722c85",
        X"a2bfe8a1", X"a81a664b", X"c24b8b70", X"c76c51a3", X"d192e819", X"d6990624", X"f40e3585", X"106aa070",
        X"19a4c116", X"1e376c08", X"2748774c", X"34b0bcb5", X"391c0cb3", X"4ed8aa4a", X"5b9cca4f", X"682e6ff3",
        X"748f82ee", X"78a5636f", X"84c87814", X"8cc70208", X"90befffa", X"a4506ceb", X"bef9a3f7", X"c67178f2");
    -- functions used in SHA
    function ch (x, y, z : std_logic_vector(0 to 32-1)) return std_logic_vector;
    function maj (x, y, z : std_logic_vector(0 to 32-1)) return std_logic_vector;
    function ss0 (x : std_logic_vector(0 to 32-1)) return std_logic_vector;
    function ss1 (x : std_logic_vector(0 to 32-1)) return std_logic_vector;
    function s0 (x : std_logic_vector(0 to 32-1)) return std_logic_vector;
    function s1 (x : std_logic_vector(0 to 32-1)) return std_logic_vector;
    -- a single round of SHA with pervious hash, round constant (kt), and message schedule (wt) input
    function sha_round (hash_prev : std_logic_vector(0 to 256 - 1);
                            kt, wt : std_logic_vector(0 to 32-1)) return std_logic_vector;
end package sha2_256_pkg;

package body sha2_256_pkg is
    function ch (x, y, z : std_logic_vector(0 to 32-1)) return std_logic_vector is
    begin
        return (x and y) xor ((not x) and z);
    end function ch;
    function maj (x, y, z : std_logic_vector(0 to 32-1)) return std_logic_vector is
    begin
        return (x and y) xor (x and z) xor (y and z);
    end function maj;
    function ss0 (x : std_logic_vector(0 to 32-1)) return std_logic_vector is
    begin
        return (x ror 2) xor (x ror 13) xor (x ror 22);
    end function ss0;
    function ss1 (x : std_logic_vector(0 to 32-1)) return std_logic_vector is
    begin
        return (x ror 6) xor (x ror 11) xor (x ror 25);
    end function ss1;
    function s0 (x : std_logic_vector(0 to 32-1)) return std_logic_vector is
    begin
        return (x ror 7) xor (x ror 18) xor (x srl 3);
    end function s0;
    function s1 (x : std_logic_vector(0 to 32-1)) return std_logic_vector is
    begin
        return (x ror 17) xor (x ror 19) xor (x srl 10);
    end function s1;
    function sha_round (hash_prev : std_logic_vector(0 to 256 - 1);
                        kt, wt : std_logic_vector(0 to 32-1)) return std_logic_vector is
        -- internal variables, 8+2 32-bit words and intermediates
        variable a, b, c, d, e, f, g, h, t1, t2 : std_logic_vector(0 to 32-1);
    begin
        a := hash_prev(0 to 32 - 1);
        b := hash_prev(32 to 64 - 1);
        c := hash_prev(64 to 96 - 1);
        d := hash_prev(96 to 128 - 1);
        e := hash_prev(128 to 160 - 1);
        f := hash_prev(160 to 192 - 1);
        g := hash_prev(192 to 224 - 1);
        h := hash_prev(224 to 256 - 1);
        t1 := std_logic_vector(unsigned(h) + unsigned(ss1(e)) + unsigned(ch(e,f,g)) + unsigned(kt) + unsigned(wt));
        t2 := std_logic_vector(unsigned(ss0(a)) + unsigned(maj(a,b,c)));
        return  (std_logic_vector(unsigned(t1) + unsigned(t2)))
                &(a)
                &(b)
                &(c)
                &(std_logic_vector(unsigned(d) + unsigned(t1)))
                &(e)
                &(f)
                &(g);
    end function sha_round;
end package body sha2_256_pkg;
