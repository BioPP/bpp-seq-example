DIRS = ExAlphabet ExSequence ExContainer
MAKE = make

all:
	-for d in $(DIRS); do (echo Building in $$d; cd $$d; $(MAKE) ); done

clean:
	-for d in $(DIRS); do (echo Cleaning in $$d; cd $$d; $(MAKE) clean ); done
