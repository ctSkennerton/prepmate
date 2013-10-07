CXX=g++
LIBS=-lz
CXXFLAGS=-Wall -O2 -DSEQAN_HAS_ZLIB
EXECUTABLE=prepmate
SEQAN=third_party/seqan-1.4.1/include

ifdef BZIP2
	override CXXFLAGS := $(CXXFLAGS) -DSEQAN_HAS_BZIP2
	override LIBS := $(LIBS) -lbz2
endif

ifdef DEBUG
	override CXXFLAGS := $(CXXFLAGS) -ggdb
endif

all: main.o $(EXECUTABLE)

clean:
	rm main.o fileManager.o fxtract

main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) -I$(SEQAN) $< -o $@

$(EXECUTABLE): main.o
	$(CXX) $(CXXFLAGS) -I$(SEQAN) -o $(EXECUTABLE) $^ $(LIBS)

test: test.cpp
	$(CXX) $(CXXFLAGS) -I$(SEQAN) -o $@ $< $(LIBS)
