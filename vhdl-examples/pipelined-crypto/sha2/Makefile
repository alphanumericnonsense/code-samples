PYTHON = python3
GHDL = ghdl
FILES = sha2_256_pkg.vhdl sha2_256_blk.vhdl sha2_256_blk_tb.vhdl
NITERS = 100

sha256tb2 :
	$(PYTHON) gen_data.py $(NITERS)
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 sha2_256_blk_tb2
	$(GHDL) -r --std=08 sha2_256_blk_tb2 -gn_iters_g=$(NITERS)

sha256tb1 :
	$(PYTHON) gen_data.py $(NITERS)
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 sha2_256_blk_tb1
	$(GHDL) -r --std=08 sha2_256_blk_tb1 -ginputfilename_g=msg.txt -gdigestfilename_g=digest.txt -gn_iters_g=$(NITERS)

clean :
	$(GHDL) clean --std=08
	$(GHDL) remove --std=08
	rm -f msg.txt digest.txt
