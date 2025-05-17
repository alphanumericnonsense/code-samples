library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

entity sha2_256_blk_tb1 is
    generic(
        inputfilename_g : string := "msg.txt";
        digestfilename_g : string := "digest.txt";
        n_iters_g : integer := 1000
    );
end entity sha2_256_blk_tb1;

architecture test of sha2_256_blk_tb1 is
    signal clk : std_logic := '0';
    signal msg_s : std_logic_vector(0 to 512-1) := (others => '0');
    signal digest_s : std_logic_vector(0 to 256-1) := (others => '0');
    signal arstn_s : std_logic := '0';
    function pad_f(m : std_logic_vector) return std_logic_vector is
        variable p : std_logic_vector(0 to 512-1) := (others => '0');
    begin
        p(0 to m'length - 1) := m;
        p(m'length) := '1';
        p(512 - 64 to 512 - 1) := std_logic_vector(to_unsigned(m'length, 64));
        return p;
    end function pad_f;
begin
    CONTROL : process is
        file inputfile : text;
        variable inputline_v : line;
        file digestfile : text;
        variable digestline_v : line;
        variable msg_v : std_logic_vector(0 to 256-1) := (others =>'0');
        variable digest_v : std_logic_vector(0 to 256-1) := (others =>'0');

        constant offset : natural := 68;
        type log_t is array(0 to offset-1) of std_logic_vector(0 to 256-1);
        variable log_s : log_t := (others => (others => '0'));
    begin
        -- open files
        file_open(inputfile, inputfilename_g, read_mode);
        file_open(digestfile, digestfilename_g, read_mode);
        for i in 0 to n_iters_g-1 loop
            report "cycle " & to_string(i);
            -- read and prep data from file
            readline(inputfile, inputline_v);
            hread(inputline_v, msg_v);
            readline(digestfile, digestline_v);
            hread(digestline_v, digest_v);

            log_s(i mod offset) := digest_v;

            wait until rising_edge(clk);

            msg_s <= pad_f(msg_v);
            report HT & to_hstring(digest_s);
            report HT & to_hstring(log_s((i+1) mod offset));
        end loop;
        std.env.stop;
    end process CONTROL;

    DUT : entity work.sha2_256_blk
        port map(
            clk_i => clk,
            msg_i => msg_s,
            digest_o => digest_s,
            arstn_i => arstn_s
        );

    CLOCK : process is

    begin
        clk <= not clk;
        wait for 1 ns;
    end process CLOCK;

    arstn_s <= '1';
end architecture test;

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

entity sha2_256_blk_tb2 is
    generic(
        n_iters_g : integer := 1000
    );
end entity sha2_256_blk_tb2;

architecture test of sha2_256_blk_tb2 is
    signal clk : std_logic := '0';
    signal msg_s : std_logic_vector(0 to 512-1) := (others => '0');
    signal digest_s : std_logic_vector(0 to 256-1) := (others => '0');
    signal arstn_s : std_logic := '0';
    function pad_f(m : std_logic_vector) return std_logic_vector is
        variable p : std_logic_vector(0 to 512-1) := (others => '0');
    begin
        p(0 to m'length - 1) := m;
        p(m'length) := '1';
        p(512 - 64 to 512 - 1) := std_logic_vector(to_unsigned(m'length, 64));
        return p;
    end function pad_f;
begin
    CONTROL : process is
    begin
        for i in 0 to n_iters_g-1 loop
            report "cycle " & to_string(i);
            wait until rising_edge(clk);
            -- msg_s <= X"61626380000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000018";
            -- report to_hstring(digest_s) & " " & "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
            msg_s <= pad_f(X"0000000000000000000000000000000000000000000000000000000000000000");
            report to_hstring(digest_s) & " " & "66687aadf862bd776c8fc18b8e9f8e20089714856ee233b3902a591d0d5f2925";
        end loop;
        std.env.stop;
    end process CONTROL;

    DUT : entity work.sha2_256_blk
        port map(
            clk_i => clk,
            msg_i => msg_s,
            digest_o => digest_s,
            arstn_i => arstn_s
        );

    CLOCK : process is

    begin
        clk <= not clk;
        wait for 1 ns;
    end process CLOCK;

    arstn_s <= '1';
end architecture test;
