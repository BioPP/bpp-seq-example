DIRS = ExAlphabet ExSequence ExContainer
MAKE = make
BIOPP_PATH = /tmp/bpp-crash-test/

all:
	-for d in $(DIRS); do (echo Building in $$d; cd $$d; $(MAKE) BIOPP_PATH=$(BIOPP_PATH)); done

clean:
	-for d in $(DIRS); do (echo Cleaning in $$d; cd $$d; $(MAKE) clean BIOPP_PATH=$(BIOPP_PATH)); done

apidoc:
	doxygen Doxyfile
