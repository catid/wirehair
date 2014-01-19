# Change your compiler settings here

# Clang seems to produce faster code
#CCPP = g++
#CC = gcc
#OPTFLAGS = -O3 -fomit-frame-pointer -funroll-loops
CCPP = clang++ -m64
CC = clang -m64
OPTFLAGS = -O4
DBGFLAGS = -g -O0 -DDEBUG
CFLAGS = -Wall -fstrict-aliasing -I./libcat -I./include
OPTLIBNAME = bin/libwirehair.a
DBGLIBNAME = bin/libwirehair_debug.a


# Object files

library_o = wirehair.o wirehair_codec.o memxor.o EndianNeutral.o

test_o = wirehair_test.o Clock.o


# Release target (default)

release : CFLAGS += $(OPTFLAGS)
release : LIBNAME = $(OPTLIBNAME)
release : library


# Debug target

debug : CFLAGS += $(DBGFLAGS)
debug : LIBNAME = $(DBGLIBNAME)
debug : library


# Library target

library : CFLAGS += $(OPTFLAGS)
library : $(library_o)
	ar rcs $(LIBNAME) $(library_o)


# test executable

test : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
test : clean $(test_o) library
	$(CCPP) $(test_o) -L./bin -lwirehair -o test
	./test


# test executable for mobile version

test-mobile : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
test-mobile : clean $(test_o)
	$(CCPP) $(test_o) -L./wirehair-mobile -lwirehair -o test
	./test


# LibCat objects

Clock.o : libcat/Clock.cpp
	$(CCPP) $(CFLAGS) -c libcat/Clock.cpp

EndianNeutral.o : libcat/EndianNeutral.cpp
	$(CCPP) $(CFLAGS) -c libcat/EndianNeutral.cpp

memxor.o : libcat/memxor.cpp
	$(CCPP) $(CFLAGS) -c libcat/memxor.cpp


# Library objects

wirehair.o : src/wirehair.cpp
	$(CCPP) $(CFLAGS) -c src/wirehair.cpp

wirehair_codec.o : src/wirehair_codec.cpp
	$(CCPP) $(CFLAGS) -c src/wirehair_codec.cpp


# Executable objects

wirehair_test.o : tests/wirehair_test.cpp
	$(CCPP) $(CFLAGS) -c tests/wirehair_test.cpp


# Cleanup

.PHONY : clean

clean :
	git submodule update --init
	-rm *.o test bin/*.a

