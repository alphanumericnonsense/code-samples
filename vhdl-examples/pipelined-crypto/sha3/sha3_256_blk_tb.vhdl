library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

entity sha3_256_tb1 is
    generic(
        inputfilename_g : string := "msg.txt";
        digestfilename_g : string := "digest.txt";
        n_iters_g : integer := 1000
    );
end entity sha3_256_tb1;

architecture test of sha3_256_tb1 is
    signal clk : std_logic := '0';
    signal msg_s : std_logic_vector(0 to 1600-1) := (others => '0');
    signal digest_s : std_logic_vector(0 to 1600-1) := (others => '0');
    signal arstn_s : std_logic := '0';
    function pad_f(m : std_logic_vector) return std_logic_vector is
        variable p : std_logic_vector(0 to 1600-1) := (others => '0');
    begin
        p(0 to m'length - 1) := m;
        p(m'length to m'length + 2 -1) := "01";
        p(m'length + 2) := '1';
        p(1088-1) := '1';
        return p;
    end function pad_f;
    function byte_reverse_f(a : std_logic_vector) return std_logic_vector is
        variable ar8 : std_logic_vector(0 to a'length-1) := (others => '0');
    begin
        for i in 0 to a'length/8 - 1 loop
            for j in 0 to 7 loop
                ar8(8*i + j) := a(8*i + 7 - j);
            end loop;
        end loop;
        return ar8;
    end function byte_reverse_f;
begin
    CONTROL : process is
        file inputfile : text;
        variable inputline_v : line;
        file digestfile : text;
        variable digestline_v : line;
        variable msg_v : std_logic_vector(0 to 256-1) := (others =>'0');
        variable digest_v : std_logic_vector(0 to 256-1) := (others =>'0');

        constant offset : natural := 28;
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

            msg_s <= pad_f(byte_reverse_f(msg_v));
            report to_hstring(byte_reverse_f(digest_s(0 to 256-1))) & " " & to_hstring(log_s((i+1) mod offset));
        end loop;
        std.env.stop;
    end process CONTROL;

    DUT : entity work.keccak_p_1600_24
        port map(
            clk_i => clk,
            state_i => msg_s,
            state_o => digest_s,
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

entity sha3_256_tb2 is
    generic(
        n_iters_g : integer := 1000
    );
end entity sha3_256_tb2;

architecture test of sha3_256_tb2 is
    signal clk : std_logic := '0';
    signal msg_s : std_logic_vector(0 to 1600-1) := (others => '0');
    signal digest_s : std_logic_vector(0 to 1600-1) := (others => '0');
    signal arstn_s : std_logic := '0';
    function pad_f(m : std_logic_vector) return std_logic_vector is
        variable p : std_logic_vector(0 to 1600-1) := (others => '0');
    begin
        p(0 to m'length - 1) := m;
        p(m'length to m'length + 2 -1) := "01";
        p(m'length + 2) := '1';
        p(1088-1) := '1';
        return p;
    end function pad_f;
    function byte_reverse_f(a : std_logic_vector) return std_logic_vector is
        variable ar8 : std_logic_vector(0 to a'length-1) := (others => '0');
    begin
        for i in 0 to a'length/8 - 1 loop
            for j in 0 to 7 loop
                ar8(8*i + j) := a(8*i + 7 - j);
            end loop;
        end loop;
        return ar8;
    end function byte_reverse_f;
begin
    CONTROL : process is
    begin
        for i in 0 to n_iters_g-1 loop
            report "cycle " & to_string(i);
            wait until rising_edge(clk);
            msg_s <= pad_f(byte_reverse_f("00110101"));
            report to_hstring(byte_reverse_f(digest_s(0 to 256-1))) & " " & "86bc56fc56af4c3cde021282f6b727ee9f90dd636e0b0c712a85d416c75e652d";
        end loop;
        std.env.stop;
    end process CONTROL;

    DUT : entity work.keccak_p_1600_24
        port map(
            clk_i => clk,
            state_i => msg_s,
            state_o => digest_s,
            arstn_i => arstn_s
        );

    CLOCK : process is

    begin
        clk <= not clk;
        wait for 1 ns;
    end process CLOCK;

    arstn_s <= '1';
end architecture test;
