CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c11
CPPFLAGS += -I. -Ikissfft-master
LDFLAGS ?=
LDLIBS ?= -lm

TARGET := equalizer_demo
SOURCES := main.c equalizer.c kissfft-master/kiss_fft.c kissfft-master/kiss_fftr.c
OBJECTS := $(SOURCES:.c=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

clean:
	rm -f $(TARGET) $(OBJECTS)
