CILKPP	= cilk++
LIBARG	= -O1 -g -lcilkutil
TARGET	= serialmerge
SRC	= $(addsuffix .cilk,$(TARGET))

all: $(TARGET)

$(TARGET): $(SRC)
	$(CILKPP) $(SRC) $(LIBARG) -o $@

clean:
	rm -f $(TARGET)
