# Change your compiler settings here

# Clang seems to produce faster code
#CCPP = g++
#CC = gcc
#OPTFLAGS = -O3 -fomit-frame-pointer -funroll-loops
CCPP = clang++ -m64
CC = clang -m64
OPTFLAGS = -O4
DBGFLAGS = -g -O0 -DDEBUG
CFLAGS = -Wall -fstrict-aliasing
OPTLIBNAME = libwirehair.a


# Object files

library_o = wirehair.o MemXOR.o EndianNeutral.o wirehair_codec_8.o


# Release target (default)

release : CFLAGS += $(OPTFLAGS)
release : LIBNAME = $(OPTLIBNAME)
release : library


# Library target

library : clean $(library_o)
	ar rcs $(LIBNAME) $(library_o)


# LibCat objects

Clock.o : Clock.cpp
	$(CCPP) $(CFLAGS) -c Clock.cpp

EndianNeutral.o : EndianNeutral.cpp
	$(CCPP) $(CFLAGS) -c EndianNeutral.cpp

MemXOR.o : MemXOR.cpp
	$(CCPP) $(CFLAGS) -c MemXOR.cpp


# Library objects

wirehair.o : wirehair.cpp
	$(CCPP) $(CFLAGS) -c wirehair.cpp

wirehair_codec_8.o : wirehair_codec_8.cpp
	$(CCPP) $(CFLAGS) -c wirehair_codec_8.cpp


# Cleanup

.PHONY : clean

clean :
	git submodule update --init
	-rm *.a *.o

