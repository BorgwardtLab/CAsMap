VERSION=$(shell grep "^version" setup.py | sed -e 's/[^\.0-9]//g')
PIP?=pip
PIP_INSTALL_FLAGS?="--user"
CYTHON?=cython
DISTUTILS?=python setup.py
DEBUG?=0
ifeq ($(DEBUG),0)
    # productive setup
    DISTUTILS_ACTUAL=$(DISTUTILS)
else
    # debug setup (see setup.py)
	DISTUTILS_ACTUAL=$(DISTUTILS) --debug
endif

all: wrapper


version:
	@echo $(VERSION)

wrapper:	sigpatsearch/wrapper.pyx sigpatsearch/wrapper.pxd
	$(CYTHON) --cplus sigpatsearch/wrapper.pyx
	$(DISTUTILS_ACTUAL) build_ext

clean:
	$(RM) sigpatsearch/wrapper.c
	$(RM) -r dist/
	(cd doc/; make clean)

install_requirements:
	$(PIP) install $(PIP_INSTALL_FLAGS) -r requirements.txt

install:
	$(PIP) install $(PIP_INSTALL_FLAGS) -e .

uninstall:
	-$(PIP) uninstall -y 'sigpatsearch==$(VERSION)'

test:
	$(DISTUTILS_ACTUAL) test

doc_html: install
	(cd doc/; make html)


package: all
	$(DISTUTILS) sdist

install_pkg: package
	$(PIP) install $(PIP_INSTALL_FLAGS) dist/sigpatsearch-$(VERSION).tar.gz
