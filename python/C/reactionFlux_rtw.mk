###########################################################################
## Makefile generated for MATLAB file/project 'reactionFlux'. 
## 
## Makefile     : reactionFlux_rtw.mk
## Generated on : Thu Nov 06 18:13:12 2014
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
COMPUTER                  = MACI64
MATLAB_ROOT               = /Applications/MATLAB_R2014b.app
MATLAB_BIN                = /Applications/MATLAB_R2014b.app/bin
MATLAB_ARCH_BIN           = /Applications/MATLAB_R2014b.app/bin/maci64
MASTER_ANCHOR_DIR         = 
START_DIR                 = /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/matlab
ARCH                      = maci64
RELATIVE_PATH_TO_ANCHOR   = .

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          Clang v3.1 | gmake (64-bit Mac)
# Supported Version(s):    3.1
# ToolchainInfo Version:   R2014b
# Specification Revision:  1.0
# 

#-----------
# MACROS
#-----------

ANSI_OPTS       = -fno-common -fexceptions
CPP_ANSI_OPTS   = -fno-common -fexceptions
ARCHS           = x86_64
XCODE_SDK_VER   = $(shell xcodebuild -showsdks | perl -anle 'BEGIN{@l = "";} push @l, $$F[-1] if /macosx/; END{ sort @l; $$_ = $$l[1]; s/macosx//; print $$_;}')
XCODE_SDK       = MacOSX$(XCODE_SDK_VER).sdk
XCODE_DEVEL_DIR = $(shell xcode-select -print-path)
XCODE_SDK_ROOT  = $(shell find $(XCODE_DEVEL_DIR) -name $(XCODE_SDK))

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: Clang C Compiler
CC = xcrun clang

# Linker: Clang Linker
LD = xcrun clang

# C++ Compiler: Clang C++ Compiler
CPP = xcrun clang++

# C++ Linker: Clang C++ Linker
CPP_LD = xcrun clang++

# Archiver: Clang Archiver
AR = xcrun ar

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_BIN)
MEX = $(MEX_PATH)/mex

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: GMAKE Utility
MAKE_PATH = %MATLAB%/bin/maci64
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
CFLAGS               = -c -isysroot $(XCODE_SDK_ROOT) -arch $(ARCHS) $(ANSI_OPTS) \
                       -O0
CPPFLAGS             = -c -isysroot $(XCODE_SDK_ROOT) -arch $(ARCHS) $(CPP_ANSI_OPTS) \
                       -O0
CPP_LDFLAGS          = -arch $(ARCHS) -isysroot $(XCODE_SDK_ROOT) -L"$(MATLAB_ARCH_BIN)"
CPP_SHAREDLIB_LDFLAGS  = -dynamiclib -isysroot $(XCODE_SDK_ROOT) -L"$(MATLAB_ARCH_BIN)" \
                         -Wl,$(LD_NAMESPACE) $(LD_UNDEFS)
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = -arch $(ARCHS) -isysroot $(XCODE_SDK_ROOT) -L"$(MATLAB_ARCH_BIN)"
MEX_CFLAGS           = -MATLAB_ARCH=$(ARCH) $(INCLUDES) \
                         \
                       COPTIMFLAGS="$(ANSI_OPTS)  \
                       -O0 \
                        $(DEFINES)" \
                         \
                       -silent
MEX_LDFLAGS          = LDFLAGS=='$$LDFLAGS'
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = -dynamiclib -isysroot $(XCODE_SDK_ROOT) -L"$(MATLAB_ARCH_BIN)" \
                       -Wl,$(LD_NAMESPACE) $(LD_UNDEFS)

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
MEX_EXT             = .mexmaci64
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

INCLUDES_BUILDINFO = -I$(START_DIR) -I/Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C -I/usr/include -I/usr/include/c++/4.2.1 -I$(MATLAB_ROOT)/extern/include -I$(MATLAB_ROOT)/simulink/include -I$(MATLAB_ROOT)/rtw/c/src -I$(MATLAB_ROOT)/rtw/c/src/ext_mode/common -I$(MATLAB_ROOT)/rtw/c/ert

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_STANDARD = -DMODEL=reactionFlux -DHAVESTDIO -DUSE_RTMODEL -DUNIX

DEFINES = $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/reactionFlux_rtwutil.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/reactionFlux_initialize.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/reactionFlux_terminate.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/reactionFlux.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/J.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/dFdTheta.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2fdTheta2.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/Afunc.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/dAdTheta.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/dAdPhi.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2AdPhi2.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2AdTheta2.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2AdThetadPhi.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2AdPhidTheta.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/Efunc.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/dEdTheta.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2EdTheta2.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2EdThetadPhi.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2EdPhidTheta.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/dEdPhi.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/d2EdPhi2.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/MI.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/S0.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/S20.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/SV0.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/S2V0.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/systemJacobian.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/rdivide.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/rt_nonfinite.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/rtGetNaN.c /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/rtGetInf.c

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = reactionFlux_rtwutil.o reactionFlux_initialize.o reactionFlux_terminate.o reactionFlux.o J.o dFdTheta.o d2fdTheta2.o Afunc.o dAdTheta.o dAdPhi.o d2AdPhi2.o d2AdTheta2.o d2AdThetadPhi.o d2AdPhidTheta.o Efunc.o dEdTheta.o d2EdTheta2.o d2EdThetadPhi.o d2EdPhidTheta.o dEdPhi.o d2EdPhi2.o MI.o S0.o S20.o SV0.o S2V0.o systemJacobian.o rdivide.o rt_nonfinite.o rtGetNaN.o rtGetInf.o

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


%.o : /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/%.c
	$(CC) $(CFLAGS) -o "$@" "$<"


%.o : /Users/justinfeigelman/Projects/LNA++/matlab/models/BirthDeath/C/%.cpp
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


