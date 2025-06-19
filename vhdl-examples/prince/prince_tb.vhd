library ieee;
	use ieee.std_logic_1164.all;
	use ieee.numeric_std.all;
	
	use work.prince_pkg.all;

entity prince_tb is

end entity prince_tb;

architecture test of prince_tb is
	
begin
	MAIN : process is
		variable key : std_logic_vector(128-1 downto 0);
		variable txtin : std_logic_vector(64-1 downto 0);
		variable txtout : std_logic_vector(64-1 downto 0);
		variable mprimein : std_logic_vector(64-1 downto 0) := (63 => '1', others => '0');
		variable mhatin : std_logic_vector(16-1 downto 0) := (15 => '1', others => '0');
		variable mbasein : std_logic_vector(4-1 downto 0) := (3 => '1', others => '0');
	begin
		-- report "mbase test";
		-- for i in 0 to 3 loop
			-- for j in 0 to 3 loop
				-- report to_string(i) & " " & to_string(j) & " " & to_string(mbase_f(mbasein, j));
			-- end loop;
			-- mbasein := mbasein srl 1;
		-- end loop;
		
		-- report "m0hat test";
		-- for i in 0 to 15 loop
			-- report to_string(i) & " " & to_string(m0hat_f(mhatin));
			-- mhatin := mhatin srl 1;
		-- end loop;
		
		-- report "m1hat test";
		-- for i in 0 to 15 loop
			-- report to_string(i) & " " & to_string(m1hat_f(mhatin));
			-- mhatin := mhatin srl 1;
		-- end loop;
		
		-- report "mprime test";
		-- for i in 0 to 63 loop
			-- report to_string(i) & " " & to_string(mprime_f(mprimein));
			-- mprimein := mprimein srl 1;
		-- end loop;
	
		-- report "shift rows test";
		-- mprimein := x"0123456789abcdef";
		-- report to_hstring(shift_rows_f(mprimein));
		-- report to_hstring(shift_rows_inv_f(mprimein));
		
		report "test vectors from https://eprint.iacr.org/2012/529.pdf";
		report "";
		key := (others => '0');
		txtin := (others => '0');
		report "key:    " & to_hstring(key);
		report "txtin:  " & to_hstring(txtin);
		report "txtout: " & to_hstring(prince_f(txtin, key, '0'));
		report "actual: 818665aa0d02dfda";
		report "";

		key := (others => '0');
		txtin := (others => '1');
		report "key:    " & to_hstring(key);
		report "txtin:  " & to_hstring(txtin);
		report "txtout: " & to_hstring(prince_f(txtin, key, '0'));
		report "actual: 604ae6ca03c20ad";
		report "";
		
		key(128-1 downto 64) := (others => '1');
		key(64-1 downto 0) := (others => '0');
		txtin := (others => '0');
		report "key:    " & to_hstring(key);
		report "txtin:  " & to_hstring(txtin);
		report "txtout: " & to_hstring(prince_f(txtin, key, '0'));
		report "actual: 9fb51935fc3df524";
		report "";
		
		key(128-1 downto 64) := (others => '0');
		key(64-1 downto 0) := (others => '1');
		txtin := (others => '0');
		report "key:    " & to_hstring(key);
		report "txtin:  " & to_hstring(txtin);
		report "txtout: " & to_hstring(prince_f(txtin, key, '0'));
		report "actual: 78a54cbe737bb7ef";
		report "";
		
		key(128-1 downto 64) := (others => '0');
		key(64-1 downto 0) := x"fedcba9876543210";
		txtin := x"0123456789abcdef";
		report "key:    " & to_hstring(key);
		report "txtin:  " & to_hstring(txtin);
		report "txtout: " & to_hstring(prince_f(txtin, key, '0'));
		report "actual: ae25ad3ca8fa9ccf";
		report "";
		
		report "enc/dec test";
		key := x"cafebabedeadbeeff007ba11ba5eba11";
		txtin := x"0123456789abcdef";
		report "key:    " & to_hstring(key);
		report "txtin:  " & to_hstring(txtin);
		report "encrypt";
		txtout := prince_f(txtin, key, '0');
		report "txtout: " & to_hstring(txtout);
		report "decrypt";
		txtin := prince_f(txtout, key, '1');
		report "txtin:  " & to_hstring(txtin);
		std.env.stop;
	end process MAIN;
end architecture test;