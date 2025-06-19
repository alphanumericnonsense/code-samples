SRC = prince_pkg.vhd prince_tb.vhd

syntax:
	ghdl -s --std=08 $(SRC)

clean:
	ghdl --clean --std=08
	ghdl --remove --std=08
	rm -rf wave.ghw
	
tb: syntax
	ghdl -a --std=08 $(SRC)
	ghdl -e --std=08 $(UNIT)
	ghdl -r --std=08 $(UNIT) --wave=wave.ghw