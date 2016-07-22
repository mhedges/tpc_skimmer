all:
	g++ -o skimmer -L$(HOME)/anaconda2/lib skimmer.C -O2 `root-config --cflags --libs` -lTreePlayer
