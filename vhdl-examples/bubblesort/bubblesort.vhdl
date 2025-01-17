library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package array_types is
    type slv_1d_t is array(natural range <>) of std_logic_vector;
end package array_types;

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity comp_swap_e is
    generic(
        width_g : positive
    );
    port(
        --clk_i : in std_logic := '0';
        int1_i, int2_i : in std_logic_vector(width_g - 1 downto 0) := (others => '0');
        int1_o, int2_o : out std_logic_vector(width_g - 1 downto 0) := (others => '0')
    );
end entity comp_swap_e;

architecture blah of comp_swap_e is

begin
    process(int1_i, int2_i) is

    begin
        if unsigned(int1_i) > unsigned(int2_i) then
            int2_o <= int1_i;
            int1_o <= int2_i;
        else
            int1_o <= int1_i;
            int2_o <= int2_i;
        end if;
    end process;
end architecture blah;

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.array_types.all;

entity bubblesort_e is
    generic(
        depth_g, width_g : positive
    );
    port(
        start_i, clk_i : in std_logic := '0';
        array_i : in slv_1d_t(0 to depth_g - 1)(width_g - 1 downto 0) := (others => (others => '0'));
        array_o : out slv_1d_t(0 to depth_g - 1)(width_g - 1 downto 0) := (others => (others => '0'));
        done_o : out std_logic := '0'
    );
end entity bubblesort_e;

architecture blah of bubblesort_e is
    signal temp_array_i_s : slv_1d_t(0 to depth_g - 1)(width_g - 1 downto 0) := (others => (others => '0'));
    signal temp_array_s : slv_1d_t(0 to depth_g - 1)(width_g - 1 downto 0) := (others => (others => '0'));
    signal temp_array_o_s : slv_1d_t(0 to depth_g - 1)(width_g - 1 downto 0) := (others => (others => '0'));
begin
    control : process(clk_i) is
        variable round_v : natural := 0;
    begin
        if rising_edge(clk_i) then
            if start_i then
                report "sorting started";
                round_v := 1;
                done_o <= '0';
                temp_array_i_s <= array_i;
            elsif (round_v > 0) and (round_v <= depth_g/2) then
                --report "round: "&to_string(round_v);
                round_v := round_v + 1;
                temp_array_i_s <= temp_array_o_s;
            elsif round_v > depth_g/2 then
                report "sorting finished";
                array_o <= temp_array_o_s;
                done_o <= '1';
            else
                NULL;
            end if;
        end if;
    end process control;

    gen_even : for i in 0 to depth_g/2 - 1 generate
        compswapeveni : entity work.comp_swap_e
            generic map(
                width_g => width_g
            )
            port map(
                --clk_i => clk_i,
                int1_i => temp_array_i_s(2*i),
                int2_i => temp_array_i_s(2*i+1),
                int1_o => temp_array_s(2*i),
                int2_o => temp_array_s(2*i+1)
            );
    end generate;

    gen_odd : for i in 0 to depth_g/2 - 2 generate
        compswapoddi : entity work.comp_swap_e
        generic map(
            width_g => width_g
        )
        port map(
            --clk_i => clk_i,
            int1_i => temp_array_s(2*i+1),
            int2_i => temp_array_s(2*i+2),
            int1_o => temp_array_o_s(2*i+1),
            int2_o => temp_array_o_s(2*i+2)
        );
    end generate;

    -- boundary assignments
    temp_array_o_s(0) <= temp_array_s(0);
    temp_array_o_s(depth_g-1) <= temp_array_s(depth_g-1);


end architecture blah;

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.array_types.all;
use std.textio.all;

entity bubblesort_tb_e is
    generic(
        width_g : positive := 8;
        depth_g : positive := 10;
        num_trials_g : positive := 10;
        sort_file_g : string := "sort.txt"
    );
end entity bubblesort_tb_e;

architecture blah of bubblesort_tb_e is
    signal clk_s, start_s, done_s : std_logic := '0';
    signal array_i_s, array_o_s : slv_1d_t(0 to depth_g - 1)(width_g - 1 downto 0) := (others => (others => '0'));
begin
    TB : process is
        file data_file : text;
        variable line_v : line;
        variable data_int_v : natural;
    begin
        -- open data file
        file_open(data_file, sort_file_g, read_mode);
        for t in 0 to num_trials_g - 1 loop
            report "trial " & to_string(t+1) & " of " & to_string(num_trials_g);
            -- load data
            for i in 0 to depth_g - 1 loop
                readline(data_file, line_v);
                read(line_v, data_int_v);
                array_i_s(i) <= std_logic_vector(to_unsigned(data_int_v, width_g));
            end loop;
            wait for 2 ns;

            -- show input array
            for i in 0 to depth_g - 1 loop
                report "unsorted input: "&to_string(to_integer(unsigned(array_i_s(i))));
            end loop;

            -- sort input array
            start_s <= '1';
            wait for 2 ns;
            start_s <= '0';
            wait until done_s;

            -- show sorted output array
            for i in 0 to depth_g -1 loop
                report "sorted output: "&to_string(to_integer(unsigned(array_o_s(i))));
            end loop;
        end loop;
        file_close(data_file);
        std.env.stop;
    end process TB;

    DUT : entity work.bubblesort_e
        generic map(
            width_g => width_g,
            depth_g => depth_g
        )
        port map(
            clk_i => clk_s,
            start_i => start_s,
            done_o => done_s,
            array_i => array_i_s,
            array_o => array_o_s
        );

    CLOCK : process(clk_s) is

    begin
        --report "fucking eh";
        clk_s <= not clk_s after 1 ns;
    end process CLOCK;
end architecture blah;
