# Default values
ARCH      ?= x86
DATA_TYPE ?= float
PLATFORM  ?= linux
NUM_RX    ?= 4
NUM_TX    ?= 4
NUM_SC    ?= 1024

# Valid values
VALID_ARCHS      := x86 rv rvv
VALID_DATA_TYPES := float fixed
VALID_PLATFORMS  := linux ara baremetal

# Validate inputs
ifneq ($(filter $(ARCH),$(VALID_ARCHS)),$(ARCH))
	$(error Invalid ARCH: $(ARCH). Supported: $(VALID_ARCHS))
endif
ifneq ($(filter $(DATA_TYPE),$(VALID_DATA_TYPES)),$(DATA_TYPE))
	$(error Invalid DATA_TYPE: $(DATA_TYPE). Supported: $(VALID_DATA_TYPES))
endif
ifneq ($(filter $(PLATFORM),$(VALID_PLATFORMS)),$(PLATFORM))
	$(error Invalid PLATFORM: $(PLATFORM). Supported: $(VALID_PLATFORMS))
endif
ifneq ($(shell test $(NUM_RX) -gt 0 >/dev/null 2>&1 && echo valid),valid)
	$(error NUM_RX must be an integer > 0)
endif
ifneq ($(shell test $(NUM_TX) -gt 0 >/dev/null 2>&1 && echo valid),valid)
	$(error NUM_TX must be an integer > 0)
endif
ifneq ($(shell test $(NUM_SC) -gt 0 >/dev/null 2>&1 && echo valid),valid)
	$(error NUM_SC must be an integer > 0)
endif

# CFLAGS
CFLAGS += -DARCH_$(ARCH)
CFLAGS += -DDATA_TYPE_$(DATA_TYPE)
CFLAGS += -DPLATFORM_$(PLATFORM)
CFLAGS += -DNUM_RX=$(NUM_RX)
CFLAGS += -DNUM_TX=$(NUM_TX)
CFLAGS += -DNUM_SC=$(NUM_SC)

# Compiler selection
ifeq ($(ARCH),x86)
	CC := gcc
else
	CC := riscv64-unknown-elf-gcc
endif

# Output file
OUTPUT := build/mmse_$(ARCH)_$(DATA_TYPE)_$(PLATFORM)_$(NUM_RX)x$(NUM_TX)x$(NUM_SC).elf

# Source files
SRCS := main.c $(wildcard src/*.c)

# Phony targets
.PHONY: all help gen_data

# Default target
all: gen_data $(OUTPUT)

# Run data generation
gen_data:
	python script/gen_data.py $(NUM_TX) $(NUM_RX) $(NUM_SC)

# Compile
$(OUTPUT): $(SRCS)
	mkdir -p build
	$(CC) $(CFLAGS) $^ -o $@

# Help target
help:
	@echo "Usage:"
	@echo "  make [ARCH=<arch>] [DATA_TYPE=<type>] [PLATFORM=<platform>] [NUM_RX=<num_rx>] [NUM_TX=<num_tx>] [NUM_SC=<num_sc>]"
	@echo ""
	@echo "Supported ARCH values:"
	@echo "  - x86 (default)"
	@echo "  - rv"
	@echo "  - rvv"
	@echo ""
	@echo "Supported DATA_TYPE values:"
	@echo "  - float (default)"
	@echo "  - fixed"
	@echo ""
	@echo "Supported PLATFORM values:"
	@echo "  - linux (default)"
	@echo "  - ara"
	@echo "  - baremetal"
	@echo ""
	@echo "Supported NUM_RX values: integers > 0 (default = 4)"
	@echo "Supported NUM_TX values: integers > 0 (default = 4)"
	@echo "Supported NUM_SC values: integers > 0 (default = 1024)"
