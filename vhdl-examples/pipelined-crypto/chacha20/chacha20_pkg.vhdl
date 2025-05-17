library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
-- https://datatracker.ietf.org/doc/html/rfc8439
package chacha20_pkg is
    type chacha20_state_t is array(0 to 3, 0 to 3) of std_logic_vector(0 to 32-1);
    type chacha20_quad_t is array(0 to 3) of std_logic_vector(0 to 32-1);
    constant r1_c : chacha20_quad_t := (X"61707865",X"3320646e",X"79622d32",X"6b206574");

    function quarter_round_f(q_i : chacha20_quad_t) return chacha20_quad_t;
    function col_round_f(s : chacha20_state_t) return chacha20_state_t;
    function diag_round_f(s : chacha20_state_t) return chacha20_state_t;
    function add_states_f(state1, state2 : chacha20_state_t) return chacha20_state_t;
    function serialize_state_f(state_i : chacha20_state_t) return std_logic_vector;
    function biglittleswap_f(wi : std_logic_vector(0 to 32-1)) return std_logic_vector;
end package chacha20_pkg;

package body chacha20_pkg is
    function quarter_round_f(q_i : chacha20_quad_t) return chacha20_quad_t is
        variable a,b,c,d : std_logic_vector(0 to 32-1) := (others => '0');
        variable q_o : chacha20_quad_t := (others => (others => '0'));
    begin
        a := q_i(0);
        b := q_i(1);
        c := q_i(2);
        d := q_i(3);
        --
        a := std_logic_vector(unsigned(a) + unsigned(b));
        d := d XOR a;
        d := d rol 16;
        --
        c := std_logic_vector(unsigned(c) + unsigned(d));
        b := b XOR c;
        b := b rol 12;
         --
        a := std_logic_vector(unsigned(a) + unsigned(b));
        d := d XOR a;
        d := d rol 8;
        --
        c := std_logic_vector(unsigned(c) + unsigned(d));
        b := b XOR c;
        b := b rol 7;
        --
        q_o := (a,b,c,d);
        return q_o;
    end function quarter_round_f;

    function col_round_f(s : chacha20_state_t) return chacha20_state_t is
        variable col_v, new_col_v : chacha20_quad_t := (others => (others => '0'));
        variable t : chacha20_state_t := (others => (others => (others => '0')));
    begin
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                col_v(j) := s(j,i);
            end loop;
            new_col_v := quarter_round_f(col_v);
            for j in 0 to 3 loop
                t(j,i) := new_col_v(j);
            end loop;
        end loop;
        return t;
    end function col_round_f;

    function diag_round_f(s : chacha20_state_t) return chacha20_state_t is
        variable col_v, new_col_v : chacha20_quad_t := (others => (others => '0'));
        variable t : chacha20_state_t := (others => (others => (others => '0')));
    begin
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                col_v(j) := s(j, (j+i) mod 4);
            end loop;
            new_col_v := quarter_round_f(col_v);
            for j in 0 to 3 loop
                t(j,(j+i) mod 4) := new_col_v(j);
            end loop;
        end loop;
        return t;
    end function diag_round_f;

    function add_states_f(state1, state2 : chacha20_state_t) return chacha20_state_t is
        variable state_o : chacha20_state_t := (others => (others => (others =>'0')));
    begin
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                state_o(i,j) := std_logic_vector(unsigned(state1(i,j)) + unsigned(state2(i,j)));
            end loop;
        end loop;
        return state_o;
    end function add_states_f;

    function serialize_state_f(state_i : chacha20_state_t) return std_logic_vector is
        variable ss : std_logic_vector(0 to 512-1) := (others => '0');
    begin
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                ss(32*(4*i+j) to 32*(4*i+j+1)-1) := biglittleswap_f(state_i(i,j));
            end loop;
        end loop;
        return ss;
    end function serialize_state_f;

    function biglittleswap_f(wi : std_logic_vector(0 to 32-1)) return std_logic_vector is
        variable wo : std_logic_vector(0 to 32-1) := (others => '0');
    begin
        for i in 0 to 3 loop
            wo(8*i to 8*(i+1)-1) := wi(8*(3-i) to 8*(3-i+1)-1);
        end loop;
        return wo;
    end function biglittleswap_f;
end package body chacha20_pkg;
