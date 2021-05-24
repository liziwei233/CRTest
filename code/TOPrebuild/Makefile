TOP=$(shell pwd)
SRC=$(TOP)/src
INCLUDE=$(TOP)/include
LIB=$(TOP)/lib

target: DircRec DircSep
OBJ=$(patsubst %.C,$(LIB)/lib%.so,$(notdir $(wildcard $(SRC)/*.C)))
libso=$(patsubst %.C,-l%,$(notdir $(wildcard $(SRC)/*.C)))

DircRec: $(LIB)/DircRec.o $(OBJ)
	@echo "Building "$(notdir $@)
	@g++ `root-config --cflags --libs` $^ -o $@

DircSep: $(LIB)/DircSep.o
	@echo "Building "$(notdir $@)
	@g++ `root-config --cflags --libs` $^ -o $@

$(LIB)/lib%.so: $(SRC)/%.C
	@echo "Building "$(notdir $@)
	@mkdir -p lib; g++ -I$(INCLUDE) `root-config --cflags --libs` $^ -fPIC -shared -o $@

$(LIB)/%.o: %.C
	@echo "Building "$(notdir $@)
	@mkdir -p lib; g++ -I$(INCLUDE) `root-config --cflags --libs` -c $^ -o $@

make clean:
	@rm -rf DircRec DircSep lib
