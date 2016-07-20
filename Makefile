all:
	g++ -o skimmer -L$(HOME)/anaconda2/lib skimmer.C `root-config --cflags --libs` -lTreePlayer
