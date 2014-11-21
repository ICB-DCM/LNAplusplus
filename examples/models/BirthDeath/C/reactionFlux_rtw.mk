###########################################################################
## Makefile generated for MATLAB file/project 'reactionFlux'. 
## 
## Makefile     : reactionFlux_rtw.mk
## Generated on : Fri Nov 21 14:20:40 2014
## MATLAB Coder version: 2.7 (R2014b)
## 
## Build Info:
## 
## Final product: $(RELATIVE_PATH_TO_ANCHOR)/reactionFlux.a
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# COMPUTER                Computer type. See the MATLAB "computer" command.

PRODUCT_NAME              = reactionFlux
MAKEFILE                  = reactionFlux_rtw.mk
COMPUTER                  = GLNXA64
MATLAB_ROOT               = /usr/local/MATLAB/R2014b
MATLAB_BIN                = /usr/local/MATLAB/R2014b/bin
MATLAB_ARCH_BIN           = /usr/local/MATLAB/R2014b/bin/glnxa64
MASTER_ANCHOR_DIR         = 
START_DIR                 = /home/justin/workspace/LNA++/examples/models/BirthDeath/matlab
ARCH                      = glnxa64
RELATIVE_PATH_TO_ANCHOR   = .
ANSI_OPTS                 = -ansi -pedantic -Wno-long-long -fwrapv
CPP_ANSI_OPTS             = -std=c++98 -pedantic -Wno-long-long -fwrapv

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          GNU gcc/g++ v4.4.x | gmake (64-bit Linux)
# Supported Version(s):    4.4.x
# ToolchainInfo Version:   R2014b
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# ANSI_OPTS
# CPP_ANSI_OPTS

#-----------
# MACROS
#-----------

WARN_FLAGS         = -Wall -W -Wwrite-strings -Winline -Wstrict-prototypes -Wnested-externs -Wpointer-arith -Wcast-align
WARN_FLAGS_MAX     = $(WARN_FLAGS) -Wcast-qual -Wshadow
CPP_WARN_FLAGS     = -Wall -W -Wwrite-strings -Winline -Wpointer-arith -Wcast-align
CPP_WARN_FLAGS_MAX = $(CPP_WARN_FLAGS) -Wcast-qual -Wshadow

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: GNU C Compiler
CC = gcc

# Linker: GNU Linker
LD = gcc

# C++ Compiler: GNU C++ Compiler
CPP = g++

# C++ Linker: GNU C++ Linker
CPP_LD = g++

# Archiver: GNU Archiver
AR = ar

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_BIN)
MEX = $(MEX_PATH)/mex

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: GMAKE Utility
MAKE_PATH = %MATLAB%/bin/glnxa64
MAKE = $(MAKE_PATH)/gmake


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -g
C_OUTPUT_FLAG       = -o
LDDEBUG             = -g
OUTPUT_FLAG         = -o
CPPDEBUG            = -g
CPP_OUTPUT_FLAG     = -o
CPPLDDEBUG          = -g
OUTPUT_FLAG         = -o
ARDEBUG             =
STATICLIB_OUTPUT_FLAG =
MEX_DEBUG           = -g
RM                  = @rm -f
ECHO                = @echo
MV                  = @mv
RUN                 =

#----------------------------------------
# "Faster Builds" Build Configuration
#----------------------------------------

ARFLAGS              = ruvs
CFLAGS               = -c $(ANSI_OPTS) -fPIC \
                       -O0
CPPFLAGS             = -c $(CPP_ANSI_OPTS) -fPIC \
                       -O0
CPP_LDFLAGS          = -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)"
CPP_SHAREDLIB_LDFLAGS  = -shared -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -Wl,--no-undefined
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)"
MEX_CFLAGS           = -MATLAB_ARCH=$(ARCH) $(INCLUDES) \
                         \
                       COPTIMFLAGS="$(ANSI_OPTS)  \
                       -O0 \
                        $(DEFINES)" \
                         \
                       -silent
MEX_LDFLAGS          = LDFLAGS=='$$LDFLAGS'
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = -shared -Wl,-rpath,"$(MATLAB_ARCH_BIN)",-L"$(MATLAB_ARCH_BIN)" -Wl,--no-undefined

#--------------------
# File extensions
#--------------------

H_EXT               = .h
OBJ_EXT             = .o
C_EXT               = .c
EXE_EXT             =
SHAREDLIB_EXT       = .so
HPP_EXT             = .hpp
OBJ_EXT             = .o
CPP_EXT             = .cpp
EXE_EXT             =
SHAREDLIB_EXT       = .so
STATICLIB_EXT       = .a
MEX_EXT             = .mexa64
MAKE_EXT            = .mk


###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = $(RELATIVE_PATH_TO_ANCHOR)/reactionFlux.a
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = -I$(START_DIR) -I/home/justin/workspace/LNA++/examples/models/BirthDeath/C -I/usr/include -I/usr/include/c++/4.2.1 -I$(MATLAB_ROOT)/extern/include -I$(MATLAB_ROOT)/simulink/include -I$(MATLAB_ROOT)/rtw/c/src -I$(MATLAB_ROOT)/rtw/c/src/ext_mode/common -I$(MATLAB_ROOT)/rtw/c/ert

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_STANDARD = -DMODEL=reactionFlux -DHAVESTDIO -DUSE_RTMODEL -DUNIX

DEFINES = $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = /home/justin/workspace/LNA++/examples/models/BirthDeath/C/reactionFlux_rtwutil.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/reactionFlux_initialize.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/reactionFlux_terminate.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/reactionFlux.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/J.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/dFdTheta.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2fdTheta2.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/Afunc.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/dAdTheta.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/dAdPhi.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2AdPhi2.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2AdTheta2.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2AdThetadPhi.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2AdPhidTheta.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/Efunc.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/dEdTheta.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2EdTheta2.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2EdThetadPhi.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2EdPhidTheta.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/dEdPhi.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/d2EdPhi2.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/MI.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/S0.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/S20.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/SV0.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/S2V0.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/Y0.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/V0.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/systemJacobian.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/rdivide.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/rt_nonfinite.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/rtGetNaN.c /home/justin/workspace/LNA++/examples/models/BirthDeath/C/rtGetInf.c

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = reactionFlux_rtwutil.o reactionFlux_initialize.o reactionFlux_terminate.o reactionFlux.o J.o dFdTheta.o d2fdTheta2.o Afunc.o dAdTheta.o dAdPhi.o d2AdPhi2.o d2AdTheta2.o d2AdThetadPhi.o d2AdPhidTheta.o Efunc.o dEdTheta.o d2EdTheta2.o d2EdThetadPhi.o d2EdPhidTheta.o dEdPhi.o d2EdPhi2.o MI.o S0.o S20.o SV0.o S2V0.o Y0.o V0.o systemJacobian.o rdivide.o rt_nonfinite.o rtGetNaN.o rtGetInf.o

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

LIBS = 

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = -lm

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CFLAGS += $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CPPFLAGS += $(CPPFLAGS_BASIC)

###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute


all : build
	@echo "### Successfully generated all binary outputs."


build : prebuild $(PRODUCT)


prebuild : 


download : build


execute : download


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	@echo "### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS)  $(PRODUCT) $(OBJS)
	@echo "### Created: $(PRODUCT)"


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

%.o : %.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : %.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(MATLAB_ROOT)/rtw/c/src/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(MATLAB_ROOT)/rtw/c/src/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : $(START_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


%.o : /home/justin/workspace/LNA++/examples/models/BirthDeath/C/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : /home/justin/workspace/LNA++/examples/models/BirthDeath/C/%.cpp
	$(CPP) $(CPPFLAGS) -o "$@" "$<"


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : $(MAKEFILE) rtw_proj.tmw


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@echo "### PRODUCT = $(PRODUCT)"
	@echo "### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@echo "### BUILD_TYPE = $(BUILD_TYPE)"
	@echo "### INCLUDES = $(INCLUDES)"
	@echo "### DEFINES = $(DEFINES)"
	@echo "### ALL_SRCS = $(ALL_SRCS)"
	@echo "### ALL_OBJS = $(ALL_OBJS)"
	@echo "### LIBS = $(LIBS)"
	@echo "### MODELREF_LIBS = $(MODELREF_LIBS)"
	@echo "### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@echo "### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@echo "### CFLAGS = $(CFLAGS)"
	@echo "### LDFLAGS = $(LDFLAGS)"
	@echo "### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@echo "### CPPFLAGS = $(CPPFLAGS)"
	@echo "### CPP_LDFLAGS = $(CPP_LDFLAGS)"
	@echo "### CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS)"
	@echo "### ARFLAGS = $(ARFLAGS)"
	@echo "### MEX_CFLAGS = $(MEX_CFLAGS)"
	@echo "### MEX_LDFLAGS = $(MEX_LDFLAGS)"
	@echo "### DOWNLOAD_FLAGS = $(DOWNLOAD_FLAGS)"
	@echo "### EXECUTE_FLAGS = $(EXECUTE_FLAGS)"
	@echo "### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files..."
	$(RM) $(PRODUCT)
	$(RM) $(ALL_OBJS)
	$(ECHO) "### Deleted all derived files."


