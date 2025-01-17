library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

entity aes_tb is
    generic(
        key_g : std_logic_vector(128-1 downto 0) := x"cafebabedeadbeefb01dfaceba5eba11"
    );
end entity aes_tb;

architecture test of aes_tb is
    signal clk, rstn, cipher_valid: std_logic := '0';
    signal plain, cipher : std_logic_vector(128 - 1 downto 0) := (others => '0');
begin

    TEST : process is
        variable counter : natural := 0;
        file text_v : text;
        variable line_v : line;
        variable gt : std_logic_vector(128-1 downto 0);
    begin
        rstn <= '0';
        wait until rising_edge(clk);
        rstn <= '1';
        wait until rising_edge(clk);

        report LF & "AES128 encryption, feedback starting from 0^128 under key " & LF & to_hstring(key_g);

        file_open(text_v, "gt.txt", read_mode);

        while counter < 50 loop
            readline(text_v, line_v);
            hread(line_v, gt);

            wait until cipher_valid;
            report "computed: " & HT & to_hstring(cipher);
            report "ground truth: " & HT &  to_hstring(gt);
            counter := counter + 1;
            plain <= cipher;
        end loop;

        file_close(text_v);
        std.env.stop;
    end process TEST;

    DUT : entity work.multicycle_aes128
        port map(
            clk_in => clk,
            rstn_in => rstn,
            key_in => key_g,
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
