BIOPP_PATH=$(HOME)/.local

all:
	g++ -std=c++11 -I$(BIOPP_PATH)/include -L$(BIOPP_PATH)/lib ExContainer.cpp -lbpp-seq -lbpp-core -o excontainer

clean:
	rm excontainer
