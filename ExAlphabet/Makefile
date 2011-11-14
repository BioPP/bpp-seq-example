BIOPP_PATH=$(HOME)/.local

all:
	g++ -I$(BIOPP_PATH)/include -L$(BIOPP_PATH)/lib ExAlphabet.cpp -lbpp-seq -lbpp-core -o exalphabet

clean:
	rm exalphabet