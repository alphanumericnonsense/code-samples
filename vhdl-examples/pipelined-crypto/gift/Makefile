GHDL = ghdl
FILES = gift128_pkg.vhdl gift128_enc.vhdl gift128_enc_tb.vhdl
NITERS = 100

gift128enc:
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 gift128_enc_tb
	$(GHDL) -r --std=08 gift128_enc_tb

clean :
	$(GHDL) clean --std=08
	$(GHDL) remove --std=08
