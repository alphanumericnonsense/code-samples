PYTHON = python3
GHDL = ghdl
FILES = aes_pkg.vhdl aes128_encrypt.vhdl aes128_encrypt_tb.vhdl
NITERS = 50

aestb :
	$(PYTHON) testvectors.py $(NITERS)
	$(GHDL) syntax --std=08 $(FILES)
	$(GHDL) -a --std=08 $(FILES)
	$(GHDL) -e --std=08 aes128_encrypt_tb
	$(GHDL) -r --std=08 aes128_encrypt_tb -ginputfilename_g=key-plain.txt -gcipherfilename_g=cipher.txt -gn_iters_g=$(NITERS)

clean :
	$(GHDL) clean --std=08
	$(GHDL) remove --std=08
	rm -rf key-plain.txt cipher.txt
