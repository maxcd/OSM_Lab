# makefile: makes thee pi program

# set the target name
all: BS_asian.exe BS_asian_rand.exe

# compile the pi.cpp file
BS_asian.exe : BS_asian.cpp
	g++ BS_asian.cpp -fopenmp -o BS_asian.exe

BS_asian_rand.exe : BS_asian_rand.cpp
	g++ BS_asian_rand.cpp -fopenmp -o BS_asian_rand.exe

# delete the exec file after execution?
clean :
	rm *.exe

