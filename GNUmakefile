CC = gcc
CFLAGS = -Iinc -fno-builtin -Wall -Wextra
OBJDIR = obj
BUILDDIR = build

SRC = $(wildcard src/*.c)
OBJ = $(SRC:src/%.c=$(OBJDIR)/%.o)
OBJ += $(OBJDIR)/main.o

TARGET = $(BUILDDIR)/mmse

# Debug target
dbg: CFLAGS += -g
dbg: $(TARGET)

all: $(TARGET)

# Compile the elf
$(TARGET): $(OBJ)
	@mkdir -p $(BUILDDIR)
	$(CC) -o $@ $^

# Compile the src/*.c object files
$(OBJDIR)/%.o: src/%.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile the main.c object files
$(OBJDIR)/main.o: main.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/* $(BUILDDIR)/*

.PHONY: all dbg clean

