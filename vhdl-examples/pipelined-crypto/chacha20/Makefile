PYTHON = python3
GHDL = ghdl
FILES = chacha20_pkg.vhdl chacha20_block.vhdl chacha20_blk_tb.vhdl
NITERS = 30

chacha20blktb2 :
	$(PYTHON) chacha20blktb2.py $(NITERS)
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 chacha20_blk_tb2
	$(GHDL) -r --std=08 chacha20_blk_tb2 -ginputfilename_g=key-nonce.txt -goutputfilename_g=state.txt -gn_iters_g=$(NITERS)

chacha20blktb1 :
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 chacha20_blk_tb1
	$(GHDL) -r --std=08 chacha20_blk_tb1

clean :
	$(GHDL) clean --std=08
	$(GHDL) remove --std=08
	rm key-nonce.txt state.txt
