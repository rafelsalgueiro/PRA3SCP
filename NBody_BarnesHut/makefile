#make file for my serial version of the project
BUILD := gui

CC=gcc

# add D_GLFW_SUPPORT define for GUI based on GLFW

cxxflags.gui = -pg -O2 -c -Wall -D D_GLFW_SUPPORT
cxxflags.no_gui = -pg -O2 -c -Wall
CFLAGS = ${cxxflags.${BUILD}}
lddflags.gui = -lm -lGL -lGLU -lglfw3 -lX11 -lXxf86vm -lpthread -ldl -lXrandr -lXi -lXinerama -lXcursor
lddflags.no_gui = -lm -lpthread -ldl
LDFLAGS = ${lddflags.${BUILD}}
SOURCES=nbfast.c
OBJECTS=nbfast.o
exec.gui = NBody_GUI
exec.no_gui = NBody
EXECUTABLE = ${exec.${BUILD}}

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -pg $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE);
