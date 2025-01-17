PYTHON = python3
GHDL = ghdl
FILES = sha3_pkg.vhdl keccak_p_1600_24.vhdl sha3_256_blk_tb.vhdl
NITERS = 100

sha3_256tb2 :
	$(PYTHON) gen_data.py $(NITERS) 256
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 sha3_256_tb2
	$(GHDL) -r --std=08 sha3_256_tb2 -gn_iters_g=$(NITERS)

sha3_256tb1 :
	$(PYTHON) gen_data.py $(NITERS) 256
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 sha3_256_tb1
	$(GHDL) -r --std=08 sha3_256_tb1 -ginputfilename_g=msg.txt -gdigestfilename_g=digest.txt -gn_iters_g=$(NITERS)

clean :
	$(GHDL) clean --std=08
	$(GHDL) remove --std=08
	rm -f msg.txt digest.txt
