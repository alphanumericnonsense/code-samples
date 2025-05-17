library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;
use work.aes_pkg.all;

entity aes128_encrypt_tb is
    generic(
        inputfilename_g : string := "key-plain.txt";
        cipherfilename_g : string := "cipher.txt";
        n_iters_g : integer := 100
    );
end entity aes128_encrypt_tb;

architecture test of aes128_encrypt_tb is
    signal clk : std_logic := '0';
    signal key_s, plain_s, cipher_s : std_logic_vector(0 to 128-1) := (others => '0');
    signal arstn_s : std_logic := '0';
begin
    CONTROL : process is
        file inputfile : text;
        variable inputline_v : line;
        file cipherfile : text;
        variable cipherline_v : line;
        variable key_v, plain_v, cipher_v : std_logic_vector(0 to 128-1) := (others =>'0');

        constant offset : natural := 14;
        type log_t is array(0 to offset-1) of std_logic_vector(0 to 128-1);
        variable log_s : log_t := (others => (others => '0'));
    begin
        -- open files
        file_open(inputfile, inputfilename_g, read_mode);
        file_open(cipherfile, cipherfilename_g, read_mode);
        for i in 0 to n_iters_g-1 loop
            report "cycle " & to_string(i);
            -- read and prep data from file
            readline(inputfile, inputline_v);
            hread(inputline_v, key_v);
            readline(inputfile, inputline_v);
            hread(inputline_v, plain_v);
            readline(cipherfile, cipherline_v);
            hread(cipherline_v, cipher_v);

            log_s(i mod offset) := cipher_v;

            assert cipher_v = aes128_encrypt_f(plain_v, key_v) severity failure;

            wait until rising_edge(clk);

            key_s <= key_v;
            plain_s <= plain_v;
            report to_hstring(cipher_s) & " " & to_hstring(log_s((i+1) mod offset));
        end loop;
        std.env.stop;
    end process CONTROL;

    DUT : entity work.aes128_encrypt
        port map(
            clk_i => clk,
            key_i => key_s,
            plain_i => plain_s,
            cipher_o => cipher_s,
            arstn_i => arstn_s
        );

    CLOCK : process is

    begin
        clk <= not clk;
        wait for 1 ns;
    end process CLOCK;

    arstn_s <= '1';
end architecture test;
