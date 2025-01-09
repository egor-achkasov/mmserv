CC = gcc
CFLAGS = -Iinclude -fno-builtin -Wall -Wextra
OBJDIR = obj
BUILDDIR = .

SRC = $(wildcard src/*.c)
OBJ = $(SRC:src/%.c=$(OBJDIR)/%.o)
OBJ += $(OBJDIR)/main.o

TARGET = $(BUILDDIR)/mmse

all: $(TARGET)

# Debug target
dbg: CFLAGS += -g
dbg: $(TARGET)

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
	rm -rf $(OBJDIR)/* $(BUILDDIR)/mmse

.PHONY: all dbg clean

