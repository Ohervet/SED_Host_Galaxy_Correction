SHELL=/bin/sh

OperatingSystem:=$(shell uname -s)

ifeq (${OperatingSystem},Linux)
	ARCH=linux
endif

ifeq (${OperatingSystem},OSF1)
	ARCH=alphacxx6
endif

Dir=./
oDir=./
iDir=-Wall
lDir=-lm

#----------host gal-----------#

CXX=g++
Dir=./
oDir=./
iDir=-Wall
lDir=-lm

PROG1=gal 

N1=eg02

gal:
    EXOBJS1=\
	$(oDir)$(N1).o \
        
    $(PROG1):	$(EXOBJS1)
	$(CXX) -o $@ $(EXOBJS1) $(lDir)

    $(oDir)/$(N1).o: $(N1).cpp
	$(CXX) -c $< $(iDir) -o $@


