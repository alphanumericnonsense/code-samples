GHDL = ghdl
GHDLFLAGS = --std=08
OBJ = lfsr.vhdl lfsr_tb.vhdl
UNIT = lfsr_tb

analyse :
	$(GHDL) analyse $(GHDLFLAGS) $(OBJ)
	$(GHDL) --elab-run $(GHDLFLAGS) $(UNIT)

clean :
	$(GHDL) --clean $(GHDLFLAGS)
	$(GHDL) --remove $(GHDLFLAGS)
