CILKPP	= cilk++
LIBARG	= -O1 -g -lcilkutil
TARGET	= apsp
SRC	= $(addsuffix .cpp,$(TARGET))

all: $(TARGET)

$(TARGET): $(SRC)
	$(CILKPP) $(SRC) $(LIBARG) -o $@

clean:
	rm -f $(TARGET)