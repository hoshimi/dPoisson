OBJS = iteration.o myUtils.o
DD = dmd
DFLGS = -wi

.SUFFIXES: .d .o
sor: sor.d $(OBJS)
	$(DD) $(DFLGS) sor.d $(OBJS)

.d.o:
	$(DD) $(DFLGS) -c $<
