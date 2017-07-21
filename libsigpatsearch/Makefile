include ../Makefile.shared

SRCS=$(shell find . -name "*.cpp" | sed -e "s/^.\///")
OBJS=$(subst .cpp,.o,$(SRCS))

all: libsigpatsearch.a libsigpatsearch.$(SHAREDEXT)

clean:
	-$(RM) $(OBJS) libsigpatsearch.a libsigpatsearch.so
	-@echo ' '

# Beware: gcc+5 (5.4.0) on Ubuntu (16.04) needs to have -fPIC as a first flag
#         when building shared library (-shared), otherwise it will fail

libsigpatsearch.a: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC Archiver'
	$(AR) -r  "libsigpatsearch.a" $(OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

libsigpatsearch.$(SHAREDEXT): $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC linker'
	$(CPP) -fPIC $(CFLAGS) $(CPPFLAGS) -o libsigpatsearch.$(SHAREDEXT) -shared $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Beware: makedepend does not work w/ C++ libs, hence the sed post-processing:
#         1. remove all /usr/... dependences, then
#         2. remove remaining empty depencies lines (w/ no deps after colon)
depend:
	makedepend -- $(CPPFLAGS) -- $(SRCS) -f- | sed -e "s/ \/usr\/\([^ ]\)*//g" | sed -e  "/:\( \)*$$/d" > Makefile.depend

%.o: %.cpp
	$(CPP) -fPIC $(CFLAGS) $(CPPFLAGS) -c $<


include Makefile.depend
