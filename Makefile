include Makefile.shared
# call w/ DEBUG=1 to enable extra printing and debugging symbols, as well as to
# disable code optimisation

all:: libs exes

libs:
	(cd libsigpatsearch; make depend all)

exes:
	(cd executables; make all)

clean:
	(cd libsigpatsearch; make clean)
	(cd executables; make clean)

doc:
	(cd doc; make)

clean_doc:
	(cd doc; make clean)
