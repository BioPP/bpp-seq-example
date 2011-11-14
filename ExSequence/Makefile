BIOPP_PATH=$(HOME)/.local

all:
	g++ -I$(BIOPP_PATH)/include -L$(BIOPP_PATH)/lib ExSequence.cpp -lbpp-seq -lbpp-core -o exsequence

clean:
	rm exsequence
