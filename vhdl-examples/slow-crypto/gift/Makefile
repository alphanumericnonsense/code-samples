GHDL = ghdl
GHDLFLAGS = --std=08
OBJ = gift128_pkg.vhdl multicycle_gift128.vhdl gift_tb.vhdl
UNIT = gift_tb

analyse :
	$(GHDL) analyse $(GHDLFLAGS) $(OBJ)
	$(GHDL) --elab-run $(GHDLFLAGS) $(UNIT)

clean :
	$(GHDL) --clean $(GHDLFLAGS)
	$(GHDL) --remove $(GHDLFLAGS)

