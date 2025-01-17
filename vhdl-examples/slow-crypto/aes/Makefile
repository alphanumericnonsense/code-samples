GHDL = ghdl
GHDLFLAGS = --std=08
OBJ = aes128_enc_pkg.vhdl multicycle_aes128.vhdl aes_tb.vhdl
UNIT = aes_tb
PYTHON = python3
SCRIPT = gt.py

analyse :
	$(PYTHON) $(SCRIPT)
	$(GHDL) analyse $(GHDLFLAGS) $(OBJ)
	$(GHDL) --elab-run $(GHDLFLAGS) $(UNIT)

clean :
	$(GHDL) --clean $(GHDLFLAGS)
	$(GHDL) --remove $(GHDLFLAGS)
	rm -f gt.txt
