library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.sha2_256_pkg.all;

entity sha2_256_blk is
    port(
        clk_i : in std_logic;
        arstn_i : in std_logic;
        msg_i : in std_logic_vector(0 to 512-1); -- already padded!
        digest_o : out std_logic_vector(0 to 256-1)
    );
end entity sha2_256_blk;

architecture pipe of sha2_256_blk is
    type pipe_state_t is array(0 to 65-1) of std_logic_vector(0 to 256-1);
    signal pipe_data : pipe_state_t := (others => (others => '0'));
    type msg_sch_t is array(0 to 64-1) of std_logic_vector(0 to 16*32-1);
    signal msg_sch_16 : msg_sch_t := (others => (others => '0'));
begin
    PIPE : process(clk_i, arstn_i) is
        variable digest_v : std_logic_vector(0 to 256-1) := (others => '0');
    begin
        if NOT arstn_i then
            pipe_data <= (others => (others => '0'));
            msg_sch_16 <= (others => (others => '0'));
            digest_o <= (others => '0');
        elsif rising_edge(clk_i) then

            -- round functions
            pipe_data(0) <= SHA2_256_INIT_c;
            for i in 1 to 65-1 loop
                pipe_data(i) <= sha_round(pipe_data(i-1), SHA2_256_ROUND_c(i-1), msg_sch_16(i-1)(0 to 32-1));
                --report to_string(i) & " " & to_hstring(pipe_data(i));
            end loop;

            -- message schedule
            -- keep 16 words for each stage, update and rotate for each stage
            msg_sch_16(0) <= msg_i;
            for i in 1 to 64-1 loop
                msg_sch_16(i) <= msg_sch_16(i-1)(1*32 to 16*32-1) &
                                 std_logic_vector(unsigned(s1(msg_sch_16(i-1)(14*32 to 15*32-1)))
                                               + unsigned(msg_sch_16(i-1)(9*32 to 10*32-1))
                                               + unsigned(s0(msg_sch_16(i-1)(1*32 to 2*32-1)))
                                               + unsigned(msg_sch_16(i-1)(0 to 1*32-1)));
            end loop;

            -- hash update and output
            for i in 0 to 8-1 loop
                digest_v(32*i to 32*(i+1)-1) := std_logic_vector(unsigned(pipe_data(64)(32*i to 32*(i+1)-1)) + unsigned(SHA2_256_INIT_ARRAY_c(i)));
            end loop;
            digest_o <= digest_v;
            --report "in block: " & to_hstring(digest_v) & " " & to_hstring(digest_o);
        else
            NULL;
        end if;
    end process PIPE;
end architecture pipe;
