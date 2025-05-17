library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

entity chacha20_blk_tb1 is
    generic(
        key_g : std_logic_vector(256-1 downto 0) := X"0000000000000000000000000000000000000000000000000000000000000000";
        nonce_g : std_logic_vector(64-1 downto 0) := X"0000000000000000";
        ctr_g : std_logic_vector(64-1 downto 0) := X"0000000000000000";
        n_iters_g : integer := 100
    );
end entity chacha20_blk_tb1;

architecture test of chacha20_blk_tb1 is
    signal clk : std_logic := '0';
    signal state_s : std_logic_vector(0 to 512-1) := (others => '0');
    signal arstn_s : std_logic := '0';
begin
    CONTROL : process(clk) is
        variable counter : natural := 0;
    begin
        if rising_edge(clk) then
            counter := counter + 1;
            report to_hstring(state_s);
            -- should be report("76b8e0ada0f13d90405d6ae55386bd28bdd219b8a08ded1aa836efcc8b770dc7da41597c5157488d7724e03fb8d84a376a43b8f41518a11cc387b669b2ee6586");
        end if;
        if counter = n_iters_g then
            std.env.stop;
        end if;
    end process CONTROL;

    DUT : entity work.chacha20_block(fast)--work.chacha20_block(slow)
        port map(
            clk_i => clk,
            key_i => key_g,
            ctr_i => ctr_g,
            nonce_i => nonce_g,
            state_o => state_s,
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

entity chacha20_blk_tb2 is
    generic(
        inputfilename_g : string := "key-nonce.txt";
        outputfilename_g : string := "state.txt";
        n_iters_g : integer := 100;
        ctr_g : std_logic_vector(0 to 64-1) := (others => '0')
    );
end entity chacha20_blk_tb2;

architecture test of chacha20_blk_tb2 is
    signal clk : std_logic := '0';
    signal key_s : std_logic_vector(0 to 256-1) := (others => '0');
    signal nonce_s : std_logic_vector(0 to 64-1) := (others => '0');
    --signal counter_s : std_logic_vector(0 to 64-1) := (others => '0');
    signal state_s : std_logic_vector(0 to 512-1) := (others => '0');
    signal arstn_s : std_logic := '0';
begin
    CONTROL : process is
        file inputfile : text;
        variable inputline_v : line;
        file outputfile : text;
        variable outputline_v : line;
        variable key_v : std_logic_vector(256-1 downto 0) := (others =>'0');
        variable nonce_v, counter_v : std_logic_vector(64-1 downto 0) := (others =>'0');
        variable state_v : std_logic_vector(512-1 downto 0) := (others =>'0');

        constant offset : natural := 14;--24;-- 14 fast, 24 slow
        type log_t is array(0 to offset-1) of std_logic_vector(512-1 downto 0);
        variable log_s : log_t := (others => (others => '0'));
    begin
        -- open files
        file_open(inputfile, inputfilename_g, read_mode);
        file_open(outputfile, outputfilename_g, read_mode);
        for i in 0 to n_iters_g-1 loop
            report "cycle " & to_string(i);
            -- read and prep data from file
            readline(inputfile, inputline_v);
            hread(inputline_v, key_v);
            readline(inputfile, inputline_v);
            hread(inputline_v, nonce_v);
            readline(outputfile, outputline_v);
            hread(outputline_v, state_v);

            log_s(i mod offset) := state_v;

            wait until rising_edge(clk);

            key_s <= key_v;
            nonce_s <= nonce_v;
            report to_hstring(state_s);
            report to_hstring(log_s((i+1) mod offset));
        end loop;
        std.env.stop;
    end process CONTROL;

    DUT : entity work.chacha20_block(fast)--work.chacha20_block(slow)
        port map(
            clk_i => clk,
            key_i => key_s,
            ctr_i => ctr_g,
            nonce_i => nonce_s,
            state_o => state_s,
            arstn_i => arstn_s
        );

    CLOCK : process is

    begin
        clk <= not clk;
        wait for 1 ns;
    end process CLOCK;

    arstn_s <= '1';
end architecture test;
