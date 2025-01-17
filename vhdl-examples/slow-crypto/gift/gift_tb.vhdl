library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity gift_tb is

end entity gift_tb;

architecture test of gift_tb is
    signal clk, rstn, cipher_valid: std_logic := '0';
    signal key, plain, cipher : std_logic_vector(128 - 1 downto 0) := (others => '0');
begin

    TEST : process is
        variable counter : natural := 0;
        variable gt : std_logic_vector(128-1 downto 0);
    begin
        rstn <= '0';
        wait until rising_edge(clk);
        rstn <= '1';
        wait until rising_edge(clk);

        report "GIFT128 encryption";

        key <= (others => '0');
        plain <= (others => '0');
        gt := x"cd0bd738388ad3f668b15a36ceb6ff92";
        wait until cipher_valid;
        report "key: " & to_hstring(key);
        report "plain: " & to_hstring(plain);
        report "cipher computed: " & HT & to_hstring(cipher);
        report "cipher ground truth: " & HT &  to_hstring(gt);

        key <= x"fedcba9876543210fedcba9876543210";
        plain <= x"fedcba9876543210fedcba9876543210";
        gt := x"8422241a6dbf5a9346af468409ee0152";
        wait until cipher_valid;
        report "key: " & to_hstring(key);
        report "plain: " & to_hstring(plain);
        report "cipher computed: " & HT & to_hstring(cipher);
        report "cipher ground truth: " & HT &  to_hstring(gt);

        key <= x"d0f5c59a7700d3e799028fa9f90ad837";
        plain <= x"e39c141fa57dba43f08a85b6a91f86c1";
        gt := x"13ede67cbdcc3dbf400a62d6977265ea";
        wait until cipher_valid;
        report "key: " & to_hstring(key);
        report "plain: " & to_hstring(plain);
        report "cipher computed: " & HT & to_hstring(cipher);
        report "cipher ground truth: " & HT &  to_hstring(gt);
        std.env.stop;
    end process TEST;

    DUT : entity work.multicycle_gift128
        port map(
            clk_in => clk,
            rstn_in => rstn,
            key_in => key,
            plain_in => plain,
            cipher_out => cipher,
            cipher_valid_out => cipher_valid
        );

    CLOCK : process is

    begin
        clk <= NOT clk;
        wait for 1 ns;
    end process CLOCK;

end architecture test;
