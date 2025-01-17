GHDL = ghdl
GHDLFLAGS = --std=08
OBJ = bubblesort.vhdl
UNIT = bubblesort_tb_e
PYTHON = python3
SCRIPT = sort_file_gen.py

analyse :
	$(PYTHON) $(SCRIPT) 10 8 10 # depth_g, width_g, num_trials_g
	$(GHDL) analyse $(GHDLFLAGS) $(OBJ)
	$(GHDL) --elab-run $(GHDLFLAGS) $(UNIT)

clean :
	$(GHDL) --clean $(GHDLFLAGS)
	$(GHDL) --remove $(GHDLFLAGS)
	rm sort.txt
