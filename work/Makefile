
gnuplot=/usr/local/bin/gnuplot

N=36
mol=Snowman-0.637556-1.0-p2gg
mol=Snowman-0.637556-1.637556-p2gg
#mol=Snowman-0.637556-1.0-pg
empty=
space=$(empty) $(empty)
shape = $(word 1, $(subst -,$(space),$(1)) )
rad = $(word 2, $(subst -,$(space),$(1)) )
dist = $(word 3, $(subst -,$(space),$(1)) )
crys = $(word 4, $(subst -,$(space),$(1)) )


all: cell

cell: program
	./program -i out.lammpstrj
	$(gnuplot) unit_cell.plot
	$(gnuplot) -e "filename='trj_contact/0000000000'" ../gnuplot/frame.plot

program:
	$(MAKE) -C ../../
	cp ../../bin/mac/analysis .

single:
	tail -2 ../crystals/$(mol).svg | head -1 | python ../pylib/createCells.py . 0 $(call rad, $(mol)) $(call dist, $(mol)) $(call crys, $(mol)) > crys.dat
	$(gnuplot) ../crys.plot

grid:
	tail -2 ../crystals/$(mol).svg | head -1 | python ../pylib/createCells.py . $N $(call rad, $(mol)) $(call dist, $(mol)) $(call crys, $(mol)) > crys.dat
	$(gnuplot) ../crys.plot

test:
	python test.py > crys.dat
	$(gnuplot) ../crys.plot

