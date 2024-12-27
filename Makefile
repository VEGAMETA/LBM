# Компилятор и флаги
CC = gcc
CFLAGS = -O3 -std=c23 -march=native -mtune=native -flto -ffunction-sections -fdata-sections -fno-asynchronous-unwind-tables -fno-stack-protector
LDFLAGS = -Wl,--gc-sections -Wl,-s -L"D:\SDKs\raylib\lib" -lraylib -lopengl32 -lgdi32 -lwinmm

# Пути
INCDIR = ./include "D:\SDKs\raylib\include"
INCFLAGS = $(addprefix -I, $(INCDIR))

SRCDIR = .
BUILDDIR = build
TARGET = lbm_simulator

# Файлы
SRC = $(SRCDIR)/src/utils.c 
SRC += $(SRCDIR)/src/lbm_3d.c $(SRCDIR)/src/lbm_2d.c $(SRCDIR)/src/lbm_1d.c 
SRC += $(SRCDIR)/src/visualize_raylib.c 
SRC += $(SRCDIR)/main.c
OBJ = $(SRC:%.c=$(BUILDDIR)/%.o)

# Правила
all: $(BUILDDIR) $(TARGET)
	./$(TARGET)

$(BUILDDIR):
	mkdir $(BUILDDIR)
	mkdir $(BUILDDIR)\src

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean
