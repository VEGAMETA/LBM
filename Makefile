CC = gcc
CFLAGS = -O3 -march=native -mtune=native -flto -ffunction-sections -fdata-sections -fno-asynchronous-unwind-tables -fno-stack-protector
LDFLAGS = -Wl,--gc-sections -Wl,-s
SRCDIR = .
INCDIR = ./lbm
BUILDDIR = build
TARGET = lbm_simulator

SRC = $(SRCDIR)/main.c $(SRCDIR)/lbm/lbm.c
OBJ = $(SRC:%.c=$(BUILDDIR)/%.o)

all: $(BUILDDIR) $(OBJ) $(TARGET)
	rm -rf $(BUILDDIR) $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)
	mkdir -p $(BUILDDIR)/lbm

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

$(BUILDDIR)/%.o: %.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@
	
.PHONY: all
