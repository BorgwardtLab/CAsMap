VERSION=$(shell grep "^Version" DESCRIPTION | sed -e 's/[^\.0-9]//g')
MAKEVARS="$(abspath ${CURDIR}/src/Makevars)"
DEBUG?=0
SOURCES=R/wrapper.r src/sigpatsearch.cpp
INSTALL_SRC?=../sigpatsearch_$(VERSION).tar.gz
INSTALL_DEST?=$(shell Rscript -e ".libPaths()[1]" | cut -f 2 -d " ")
ifeq ($(DEBUG),0)
    # productive setup (with some extra CRAN curator's flags)
    PKG_CXXFLAGS=-mtune=native -Wall -pedantic
else
    # debug setup: stop at minor compilation errors, and add debug symbol
	# Note: "-O0 -g3" flags will have no effect as they are overriden with
	#       "-g -O2" flags by devtools induced compilation of C++ sources
	PKG_CXXFLAGS=-mtune=native -Wall -pedantic -Werror -DDEBUG
endif

define compile
	Rscript -e "Rcpp::compileAttributes()" && Rscript -e "devtools::document()"
endef

define compile_cran
	echo "PKG_CXXFLAGS=$(PKG_CXXFLAGS)" >> $(MAKEVARS)
	-$(call compile)
	sed -i -e "$$ d" $(MAKEVARS)
	@$(RM) $(MAKEVARS)-e # dirty fix for incompatibility of sed's -i switch (Mac OS vs. Unix)
endef

all: package

version:
	@echo $(VERSION)

wrapper:
	$(call compile_cran)

clean_src:
	$(RM) src/*.o src/*.so

clean_doc:
	$(RM) man/*.Rd

clean_check:
	$(RM) -r sigpatsearch.Rcheck/

clean_pkg:
	$(RM) $(INSTALL_SRC)

clean: clean_src clean_doc clean_check clean_pkg

package:
	Rscript -e "devtools::build()"

install: package
	R CMD INSTALL $(INSTALL_SRC) $(INSTALL_DEST)

install_cwd:
	R CMD INSTALL . $(INSTALL_DEST)

uninstall:
	-R CMD REMOVE sigpatsearch

check: package
	R CMD check --as-cran ../sigpatsearch_$(VERSION).tar.gz
