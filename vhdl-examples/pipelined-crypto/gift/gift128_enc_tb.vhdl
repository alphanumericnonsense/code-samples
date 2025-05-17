library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity gift128_enc_tb is

end entity gift128_enc_tb;

architecture test of gift128_enc_tb is
    signal clk, arstn : std_logic := '0';
    signal key, plain, cipher1, cipher2 : std_logic_vector(128-1 downto 0) := (others => '0');
begin

    TEST : process is

    begin
        arstn <= '0';
        wait for 2 ns;
        arstn <= '1';
        wait for 2 ns;
        report LF & "test 1, key=plain=0" & LF;
        key <= (others => '0');
        plain <= (others => '0');
        for i in 0 to 45 loop
            wait until rising_edge(clk);
            report "cycle " & to_string(i) & " slow: " & to_hstring(cipher1) & " " & " actual: cd0bd738388ad3f668b15a36ceb6ff92";
            report "cycle " & to_string(i) & " fast: " & to_hstring(cipher2) & " " & " actual: cd0bd738388ad3f668b15a36ceb6ff92";
        end loop;

        arstn <= '0';
        wait for 2 ns;
        arstn <= '1';
        wait for 2 ns;
        report LF & "test 2, key=plain=fedcba9876543210fedcba9876543210" & LF;
        key <= x"fedcba9876543210fedcba9876543210";
        plain <= x"fedcba9876543210fedcba9876543210";
        for i in 0 to 45 loop
            wait until rising_edge(clk);
            report "cycle " & to_string(i) & " slow: " & to_hstring(cipher1) & " " & " actual: 8422241a6dbf5a9346af468409ee0152";
            report "cycle " & to_string(i) & " fast: " & to_hstring(cipher2) & " " & " actual: 8422241a6dbf5a9346af468409ee0152";
        end loop;

        arstn <= '0';
        wait for 2 ns;
        arstn <= '1';
        wait for 2 ns;
        report LF & "test 3, key=d0f5c59a7700d3e799028fa9f90ad837, plain=e39c141fa57dba43f08a85b6a91f86c1" & LF;
        key <= x"d0f5c59a7700d3e799028fa9f90ad837";
        plain <= x"e39c141fa57dba43f08a85b6a91f86c1";
        for i in 0 to 45 loop
            wait until rising_edge(clk);
            report "cycle " & to_string(i) & " slow: " & to_hstring(cipher1) & " " & " actual: 13ede67cbdcc3dbf400a62d6977265ea";
            report "cycle " & to_string(i) & " fast: " & to_hstring(cipher2) & " " & " actual: 13ede67cbdcc3dbf400a62d6977265ea";
        end loop;

        std.env.stop;
    end process TEST;

    DUT1 : entity work.gift128_enc(pipe)
        port map(
            clk_i => clk,
            arstn_i => arstn,
            key_i => key,
            plain_i => plain,
            cipher_o => cipher1
        );

    DUT2 : entity work.gift128_enc(fast)
        port map(
            clk_i => clk,
            arstn_i => arstn,
            key_i => key,
            plain_i => plain,
            cipher_o => cipher2
        );

    CLOCK : process is

    begin
        clk <= not clk;
        wait for 1 ns;
    end process CLOCK;
end architecture test;
