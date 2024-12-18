# Компилятор и флаги
CC = gcc
CFLAGS = -O3 -std=c23 -march=native -mtune=native -flto -ffunction-sections -fdata-sections -fno-asynchronous-unwind-tables -fno-stack-protector
LDFLAGS = -Wl,--gc-sections -Wl,-s -L"D:\SDKs\raylib\lib" -lraylib -lopengl32 -lgdi32 -lwinmm

# Пути
INCDIR = ./lbm ./visualize "D:\SDKs\raylib\include"
INCFLAGS = $(addprefix -I, $(INCDIR))

SRCDIR = .
BUILDDIR = build
TARGET = lbm_simulator

# Файлы
SRC = $(SRCDIR)/visualize/visualize_raylib.c $(SRCDIR)/main.c $(SRCDIR)/lbm/lbm.c 
OBJ = $(SRC:%.c=$(BUILDDIR)/%.o)

# Правила
all: $(BUILDDIR) $(TARGET)
	./$(TARGET)

$(BUILDDIR):
	mkdir $(BUILDDIR)
	mkdir $(BUILDDIR)\lbm
	mkdir $(BUILDDIR)\visualize

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean
