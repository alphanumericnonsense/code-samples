--
-- PRINCE low latency (single-cycle) encryption/decryption
-- developed/used by NXP for, e.g., encrypted flash
-- https://eprint.iacr.org/2012/529.pdf
-- below might need lots of optimization to ever fit into a single cycle...
-- 

library ieee;
	use ieee.std_logic_1164.all;
	use ieee.numeric_std.all;

package prince_pkg is
	----------------
	-- constants
	----------------
	constant alpha_c : std_logic_vector(64-1 downto 0) := x"c0ac29b7c97c50dd";
	type RC_t is array(0 to 11) of std_logic_vector(64-1 downto 0);
	-- pi_hex = 3 . <64 bits> | RC(1) | ... | RC(5) | alpha
	constant RC_c : RC_t := (
		x"0000000000000000",
		x"13198a2e03707344",
		x"a4093822299f31d0",
		x"082efa98ec4e6c89",
		x"452821e638d01377",
		x"be5466cf34e90c6c",
		x"7ef84f78fd955cb1", -- RC(5) XOR alpha
		x"85840851f1ac43aa", -- RC(4) XOR alpha
		x"c882d32f25323c54", -- RC(3) XOR alpha
		x"64a51195e0e3610d", -- RC(2) XOR alpha
		x"d3b5a399ca0c2399", -- RC(1) XOR alpha
		x"c0ac29b7c97c50dd"  -- RC(0) XOR alpha
	);
	type sbox_t is array(0 to 15) of std_logic_vector(4-1 downto 0);
	constant SBOX_c : SBOX_t := (x"b", x"f", x"3", x"2", x"a", x"c", x"9", x"1", x"6", x"7", x"8", x"0", x"e", x"5", x"d", x"4");
	constant SBOX_INV_c : SBOX_t := (x"b", x"7", x"3", x"2", x"f", x"d", x"8", x"9", x"a", x"6", x"4", x"0", x"5", x"e", x"c", x"1");
	
	---------------------
	-- functions
	---------------------
	function sbox_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	function sbox_inv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	function shift_rows_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	function shift_rows_inv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	function mbase_f(nibble : std_logic_vector(4-1 downto 0); idx : natural) return std_logic_vector;
	function m0hat_f(state : std_logic_vector(16-1 downto 0)) return std_logic_vector;
	function m1hat_f(state : std_logic_vector(16-1 downto 0)) return std_logic_vector;
	function mprime_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	-- front
	function sm_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	-- back
	function minvsinv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	-- middle
	function smprimesinv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	function prince_rnd_f(
		state : std_logic_vector(64-1 downto 0);
		round : natural;
		k1 : std_logic_vector(64-1 downto 0)
	) return std_logic_vector;
	function prince_inv_rnd_f(
		state : std_logic_vector(64-1 downto 0);
		round : natural;
		k1 : std_logic_vector(64-1 downto 0)
	) return std_logic_vector;
	function prince_f(
		txt : std_logic_vector(64-1 downto 0);
		key : std_logic_vector(128-1 downto 0);
		enc_dec : std_logic
	) return std_logic_vector;
	function prince_core_f(state, k1 : std_logic_vector(64-1 downto 0)) return std_logic_vector;
	
end package prince_pkg;

package body prince_pkg is

function sbox_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	variable newstate : std_logic_vector(64-1 downto 0) := (others => '0');
	variable nibble : std_logic_vector(4-1 downto 0) := (others => '0');
begin
	for i in 0 to 15 loop
		nibble := state(4*(i+1)-1 downto 4*i);
		newstate(4*(i+1)-1 downto 4*i) := SBOX_c(to_integer(unsigned(nibble)));
	end loop;
	return newstate;
end function sbox_f;

function sbox_inv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	variable newstate : std_logic_vector(64-1 downto 0) := (others => '0');
	variable nibble : std_logic_vector(4-1 downto 0) := (others => '0');
begin
	for i in 0 to 15 loop
		nibble := state(4*(i+1)-1 downto 4*i);
		newstate(4*(i+1)-1 downto 4*i) := SBOX_INV_c(to_integer(unsigned(nibble)));
	end loop;
	return newstate;
end function sbox_inv_f;

function shift_rows_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	variable newstate : std_logic_vector(64-1 downto 0) := (others => '0');
begin
	for i in 0 to 3 loop
		for j in 0 to 3 loop
			newstate(64-(16*j+4*i) - 1 downto 64-(16*j+4*(i+1))) := 
				state(64-(16*((j+i) mod 4)+4*i) - 1 downto 64-(16*((j+i) mod 4)+4*(i+1)));
		end loop;
	end loop;
	return newstate;
end function shift_rows_f;

function shift_rows_inv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	variable newstate : std_logic_vector(64-1 downto 0) := (others => '0');
begin
	for i in 0 to 3 loop
		for j in 0 to 3 loop
			newstate(64-(16*j+4*i) - 1 downto 64-(16*j+4*(i+1))) := 
				state(64-(16*((j-i) mod 4)+4*i) - 1 downto 64-(16*((j-i) mod 4)+4*(i+1)));
		end loop;
	end loop;
	return newstate;
end function shift_rows_inv_f;

-- base 4x4 matrix multiplication, projection onto 3 of the 4 coords
function mbase_f(nibble : std_logic_vector(4-1 downto 0); idx : natural) return std_logic_vector is
	variable retval : std_logic_vector(4-1 downto 0) := (others => '0'); -- need to define output endianess...
begin
	case idx is
		when 0 =>
			retval := (3 => '0', 2 => nibble(2), 1 => nibble(1), 0 => nibble(0));
			return retval;
		when 1 =>
			retval := (3 => nibble(3), 2 => '0', 1 => nibble(1), 0 => nibble(0));
			return retval;
		when 2 =>
			retval := (3 => nibble(3), 2 => nibble(2), 1 => '0', 0 => nibble(0));
			return retval;
		when 3 =>
			retval := (3 => nibble(3), 2 => nibble(2), 1 => nibble(1), 0 => '0');
			return retval;
		when others =>
			NULL;
	end case;
end function mbase_f;

-- base0 16x16 matrix multiplication, block cyclic with 4x4 basic blocks, each block-row left-shifted by 0-3
function m0hat_f(state : std_logic_vector(16-1 downto 0)) return std_logic_vector is
	variable newstate : std_logic_vector(16-1 downto 0) := (others => '0');
begin
	newstate := (others => '0');
	for i in 0 to 3 loop
		for j in 0 to 3 loop
			newstate(16-4*i-1 downto 16-4*(i+1)) := newstate(16-4*i-1 downto 16-4*(i+1)) 
													XOR mbase_f(state(16-4*j-1 downto 16-4*(j+1)), (j+i) mod 4);
		end loop;
	end loop;
	return newstate;
end function m0hat_f;

-- base1 16x16 matrix multiplication, block cyclic with 4x4 basic blocks, each block-row left-shifted by 1-4
function m1hat_f(state : std_logic_vector(16-1 downto 0)) return std_logic_vector is
	variable newstate : std_logic_vector(16-1 downto 0) := (others => '0');
begin
	newstate := (others => '0');
	for i in 0 to 3 loop
		for j in 0 to 3 loop
			newstate(16-4*i-1 downto 16-4*(i+1)) := newstate(16-4*i-1 downto 16-4*(i+1)) 
													XOR mbase_f(state(16-4*j-1 downto 16-4*(j+1)), (j+i+1) mod 4);
		end loop;
	end loop;
	return newstate;
end function m1hat_f;

-- 64x64 matrix multiplication, block diagonal
function mprime_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is

begin
	return  m0hat_f(state(64-1 downto 48)) & 
			m1hat_f(state(48-1 downto 32)) & 
			m1hat_f(state(32-1 downto 16)) & 
			m0hat_f(state(16-1 downto 0));
end function mprime_f;

function sm_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is

begin
	return shift_rows_f(mprime_f(sbox_f(state)));
end function sm_f;

function minvsinv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	
begin
	return sbox_inv_f(mprime_f(shift_rows_inv_f(state)));
end function minvsinv_f;

function smprimesinv_f(state : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	
begin
	return sbox_inv_f(mprime_f(sbox_f(state)));
end function smprimesinv_f;

function prince_rnd_f(state : std_logic_vector(64-1 downto 0); round : natural; k1 : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	
begin
	return sm_f(state) XOR RC_c(round) XOR k1;
end function prince_rnd_f;

function prince_inv_rnd_f(state : std_logic_vector(64-1 downto 0); round : natural; k1 : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	
begin
	return minvsinv_f(state XOR RC_c(round) XOR k1);
end function prince_inv_rnd_f;

function prince_f(txt : std_logic_vector(64-1 downto 0); key : std_logic_vector(128-1 downto 0); enc_dec : std_logic) return std_logic_vector is
	variable k0, k0prime, k1 : std_logic_vector(64-1 downto 0) := (others => '0');
begin
	if enc_dec = '0' then
		k0 := key(128-1 downto 64);
		k0prime := (k0 ror 1) XOR (k0 srl 63);
		k1 := key(64-1 downto 0);
	else
		k0prime := key(128-1 downto 64);
		k0 := (k0prime ror 1) XOR (k0prime srl 63);
		k1 := key(64-1 downto 0) XOR alpha_c;
	end if;
	
	return prince_core_f(txt XOR k0, k1) XOR k0prime;
end function prince_f;

function prince_core_f(state, k1 : std_logic_vector(64-1 downto 0)) return std_logic_vector is
	variable newstate : std_logic_vector(64-1 downto 0) := (others => '0');
begin
	newstate := state XOR k1 XOR RC_c(0);
	--report "0: " & to_hstring(newstate);
	for rnd in 1 to 5 loop
		newstate := prince_rnd_f(newstate, rnd, k1);
		--report to_string(rnd) & ": " & to_hstring(newstate);
	end loop;
	
	newstate := smprimesinv_f(newstate);
	--report "mid: " & to_hstring(newstate);
	
	for rnd in 6 to 10 loop
		newstate := prince_inv_rnd_f(newstate, rnd, k1);
		--report to_string(rnd) & ": " & to_hstring(newstate);
	end loop;
	
	return newstate XOR RC_c(11) XOR k1;
end function prince_core_f;

end package body;