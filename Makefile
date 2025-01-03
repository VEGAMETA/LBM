CC = clang

CFLAGS = -std=c23 -O3 -DNDEBUG 
CFLAGS += -fopenmp -march=native -mtune=native -flto -ffunction-sections -fdata-sections -fno-asynchronous-unwind-tables -fno-stack-protector

LDFLAGS = -Wl,--gc-sections -Wl,-s
LDFLAGS += -lraylib -lopengl32 -lgdi32 -lwinmm -ljpeg -lzlib -lpng -lOpenCL

LDDIR = ./lib
LDDIR += "D:\SDKs\OpenCL\lib"
LDDIR += "D:\SDKs\raylib\lib"
LDDIR += "D:\SDKs\libjpeg\lib" # libjpeg binary dll is nesessary
LDDIR += "D:\SDKs\zlib\lib"
LDDIR += "D:\SDKs\libpng\lib"

INCDIR = ./include
INCDIR += "D:\SDKs\OpenCL\include"
INCDIR += "D:\SDKs\raylib\include"
INCDIR += "D:\SDKs\libjpeg\include"
INCDIR += "D:\SDKs\zlib\include"
INCDIR += "D:\SDKs\libpng\include"

LDFLAGS += $(addprefix -L, $(LDDIR))
INCFLAGS = $(addprefix -I, $(INCDIR))

SRCDIR = .
BUILDDIR = build
TARGET = lbm

SRC =
SRC += $(SRCDIR)/src/lbm/utils.c
SRC += $(SRCDIR)/src/lbm/lbm_3d.c $(SRCDIR)/src/lbm/lbm_2d.c $(SRCDIR)/src/lbm/lbm_1d.c 
SRC += $(SRCDIR)/src/visualization/visualize.c
SRC += $(SRCDIR)/src/misc/images_handler.c
SRC += $(SRCDIR)/main.c

OBJ = $(SRC:%.c=$(BUILDDIR)/%.o)


all: $(BUILDDIR) $(TARGET)
	./$(TARGET) "./test2.png"

$(BUILDDIR):
	mkdir $(BUILDDIR)
	mkdir $(BUILDDIR)\src
	mkdir $(BUILDDIR)\src\lbm
	mkdir $(BUILDDIR)\src\visualization
	mkdir $(BUILDDIR)\src\misc

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

$(BUILDDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR)

test:
	./$(TARGET) "./test.png"

.PHONY: all clean test
