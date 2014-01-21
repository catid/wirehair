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

library_o = wirehair.o MemXOR.o EndianNeutral.o Galois256.o

test_o = wirehair_test.o Clock.o
gf_test_o = gf_test.o Clock.o MemXOR.o


# Release target (default)

release : CFLAGS += $(OPTFLAGS)
release : LIBNAME = $(OPTLIBNAME)
release : library_o += wirehair_codec_8.o
release : wirehair_codec_8.o library

release-16 : CFLAGS += $(OPTFLAGS)
release-16 : LIBNAME = $(OPTLIBNAME)
release-16 : library_o += wirehair_codec_16.o
release-16 : wirehair_codec_16.o library


# Debug target

debug : CFLAGS += $(DBGFLAGS)
debug : LIBNAME = $(DBGLIBNAME)
debug : library_o += wirehair_codec_8.o
debug : wirehair_codec_8.o library

debug-16 : CFLAGS += $(DBGFLAGS)
debug-16 : LIBNAME = $(DBGLIBNAME)
debug-16 : library_o += wirehair_codec_16.o
debug-16 : wirehair_codec_16.o library


# Library target

library : $(library_o)
	ar rcs $(LIBNAME) $(library_o)


# test executable

test : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
test : $(test_o)
	$(CCPP) $(test_o) -L./bin -lwirehair -o test
	./test

test-debug : CFLAGS += -DUNIT_TEST $(DBGFLAGS)
test-debug : $(test_o)
	$(CCPP) $(test_o) -L./bin -lwirehair_debug -o test
	./test


# gf-test executable

gf-test : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
gf-test : clean $(gf_test_o)
	$(CCPP) $(gf_test_o) -o gf_test
	./gf_test


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

MemXOR.o : libcat/MemXOR.cpp
	$(CCPP) $(CFLAGS) -c libcat/MemXOR.cpp

Galois256.o : libcat/Galois256.cpp
	$(CCPP) $(CFLAGS) -c libcat/Galois256.cpp


# Library objects

wirehair.o : src/wirehair.cpp
	$(CCPP) $(CFLAGS) -c src/wirehair.cpp

wirehair_codec_8.o : src/wirehair_codec_8.cpp
	$(CCPP) $(CFLAGS) -c src/wirehair_codec_8.cpp

wirehair_codec_16.o : src/wirehair_codec_16.cpp
	$(CCPP) $(CFLAGS) -c src/wirehair_codec_16.cpp


# Executable objects

wirehair_test.o : tests/wirehair_test.cpp
	$(CCPP) $(CFLAGS) -c tests/wirehair_test.cpp

gf_test.o : tests/gf_test.cpp
	$(CCPP) $(CFLAGS) -c tests/gf_test.cpp


# Cleanup

.PHONY : clean

clean :
	git submodule update --init
	-rm bin/*.a test *.o

