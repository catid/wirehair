codec_objects = Wirehair.o EndianNeutral.o memxor.o
tester_objects = Tester.o Clock.o $(codec_objects)

CFLAGS = -O4 -Wall -fstrict-aliasing

tester : $(tester_objects)
	clang++ -o tester $(tester_objects)

Clock.o : Clock.cpp
	clang++ $(CFLAGS) -c Clock.cpp

Tester.o : Tester.cpp
	clang++ $(CFLAGS) -c Tester.cpp

Wirehair.o : codec_source/Wirehair.cpp
	clang++ $(CFLAGS) -c codec_source/Wirehair.cpp

EndianNeutral.o : codec_source/EndianNeutral.cpp
	clang++ $(CFLAGS) -c codec_source/EndianNeutral.cpp

memxor.o : codec_source/memxor.cpp
	clang++ $(CFLAGS) -c codec_source/memxor.cpp

.PHONY : clean

clean :
	-rm tester $(tester_objects)

