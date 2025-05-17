library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.gift128_pkg.all;

entity gift128_enc is
    port(
        clk_i : in std_logic;
        arstn_i : in std_logic;
        key_i : in std_logic_vector(128-1 downto 0);
        plain_i : in std_logic_vector(128-1 downto 0);
        cipher_o : out std_logic_vector(128-1 downto 0)
    );
end entity gift128_enc;

architecture pipe of gift128_enc is
    type buff_t is array(0 to 40) of std_logic_vector(128-1 downto 0);
    signal state_buff : buff_t := (others => (others => '0'));
    signal key_buff : buff_t := (others => (others => '0'));
begin
    MAIN : process(clk_i, arstn_i) is

    begin
        if NOT arstn_i then
            cipher_o <= (others => '0');
            state_buff <= (others => (others => '0'));
            key_buff <= (others => (others => '0'));
        elsif rising_edge(clk_i) then
            state_buff(0) <= plain_i;
            key_buff(0) <= key_i;
            for i in 1 to 40 loop
                state_buff(i) <= gift128_round_f(state_buff(i-1), key_buff(i-1)(6*16-1 downto 4*16), key_buff(i-1)(32-1 downto 0), i-1);
                key_buff(i) <= key_round_f(key_buff(i-1));
            end loop;
            cipher_o <= state_buff(40);
            -- report "state norm";
            -- for i in 0 to 40 loop
            --     report to_string(i) & " " & to_hstring(state_buff(i));
            -- end loop;
            -- report "key norm";
            -- for i in 0 to 40 loop
            --     report to_string(i) & " " & to_hstring(key_buff(i));
            -- end loop;
        else
            NULL;
        end if;
    end process MAIN;
end architecture pipe;

architecture fast of gift128_enc is
    type key_buff_t is array(0 to 20) of std_logic_vector(256-1 downto 0); -- holds two keys in half the stages
    type state_buff_t is array(0 to 20) of std_logic_vector(128-1 downto 0);
    signal state_buff : state_buff_t := (others => (others => '0'));
    signal key_buff : key_buff_t := (others => (others => '0'));
begin
    MAIN : process(clk_i, arstn_i) is

    begin
        if NOT arstn_i then
            cipher_o <= (others => '0');
            state_buff <= (others => (others => '0'));
            key_buff <= (others => (others => '0'));
        elsif rising_edge(clk_i) then
            state_buff(0) <= plain_i;
            key_buff(0) <= key_i&key_round_f(key_i);
            for i in 1 to 20 loop
                state_buff(i) <= double_gift128_round_f(state_buff(i-1), key_buff(i-1)(128+6*16-1 downto 128+4*16)&key_buff(i-1)(6*16-1 downto 4*16), key_buff(i-1)(128+32-1 downto 128+0)&key_buff(i-1)(32-1 downto 0), 2*i-2);
                key_buff(i) <= double_key_round_f(key_buff(i-1)(128-1 downto 0));
            end loop;
            cipher_o <= state_buff(20);
            -- report "state fast";
            -- for i in 0 to 20 loop
            --     report to_string(i) & " " & to_hstring(state_buff(i));
            -- end loop;
            -- report "key fast";
            -- for i in 0 to 20 loop
            --     report to_string(i) & " " & to_hstring(key_buff(i));
            -- end loop;
        else
            NULL;
        end if;
    end process MAIN;
end architecture fast;
