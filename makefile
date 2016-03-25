.DEFAULT_GOAL := LIB

CC = /opt/gcc-5.3.0/bin/g++
CFLAGS = -std=c++14 -Wall -Wextra
OBJDIR = ../
LIBDIR = ../
LIBFILE = libPotentialPack.a

CPPFILES = Potential.cpp \
	Potential_HarmonicBond.cpp \
	Potential_HarmonicAngle.cpp \
	Potential_StillingerWeber.cpp \
	Potential_LennardJones.cpp \
	Potential_Buckingham.cpp

OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(CPPFILES))


$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o
	rm -f $(OBJDIR)/console

LIB: $(OBJS)
	ar rvs $(LIBDIR)/$(LIBFILE) $(OBJS)


test: $(OBJDIR)/main.o $(OBJS)
	$(CC) $(OBJDIR)/main.o $(OBJS) -Wl,--rpath=/opt/gcc-5.3.0/lib64 -o $(OBJDIR)/console

