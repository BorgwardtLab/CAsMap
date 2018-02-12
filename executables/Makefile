include ../Makefile.shared

SRCS=$(shell find . -name "*.cpp" | sed -e "s/^.\///")
OBJS=$(subst .cpp,.o,$(SRCS))
ALL=$(subst .o,,$(OBJS))
LIBH=../libsigpatsearch
LIBS=$(addprefix $(LIBH)/,libsigpatsearch.a)

all: $(ALL)

significant_interval_search_exact: significant_interval_search_exact.o $(LIBS)
	$(CPP) -o $@ $^

significant_interval_search_chi: significant_interval_search_chi.o $(LIBS)
	$(CPP) -o $@ $^

significant_interval_search_wy_exact: significant_interval_search_wy_exact.o $(LIBS)
	$(CPP) -o $@ $^

significant_interval_search_wy_chi: significant_interval_search_wy_chi.o $(LIBS)
	$(CPP) -o $@ $^

significant_interval_search_fastcmh: significant_interval_search_fastcmh.o $(LIBS)
	$(CPP) -o $@ $^

significant_itemset_search_facs: significant_itemset_search_facs.o $(LIBS)
	$(CPP) -o $@ $^

clean:
	-$(RM) $(OBJS) $(ALL)
	-@echo ' '

%.o: %.cpp
	$(CPP) $(CFLAGS) $(CPPFLAGS) -I$(LIBH) -c $<
