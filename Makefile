# Build all source files.
#
# Targets in this file are:
# all     Compile and link all source files.
# clean   Remove all intermediate files.
# help    Display this information.
#
# Copyright (c) 2021 Rubens AMARO
# Distributed under the MIT License.

# Directories
BINARY_DIR := bin
SOURCE_DIR := src
OBJECT_DIR := obj
OBJECT_LIB_DIR := $(OBJECT_DIR)/lib
INCLUDE_DIR := include
TARGET_DIR := src
TARGET := main.cpp
LIB_DIR := lib

# C compiler (GNU/INTEL) (g++) or (icpc)
CC := g++
CXX := g++
DEBUGS := -O3
#DEBUGS := -O3 -g

# Got some preprocessor flags to pass ?
# -I is a preprocessor flag, not a compiler flag

#CXXFLAGS := $(DEBUGS) -std=c++11 -Wall -Wextra -MP -MMD
#CXXFLAGS := $(DEBUGS) -std=c++11 -Wall -Wextra -MP -MMD -fopenmp -lm -fsanitize=leak -fanalyzer
CXXFLAGS := $(DEBUGS) -std=c++11 -MP -MMD -fopenmp -lm -D_GNU_SOURCE
CPPFLAGS := -I $(INCLUDE_DIR) -Ieigen -Ijson/single_include/nlohmann -Ilibigl/include

# Got some linker flags ?
# -L is a linker flag
LDFLAGS := -L $(LIB_DIR)
#LDLIBS := -lboost_regex
LDLIBS := -lstdc++fs

MKDIR := mkdir -p
MV := mv -f
RM := rm -rf
SED := sed
TEST := test

# Creates an object directory if it does not exist.
create_binary_directory := $(shell for f in $(BINARY_DIR); do $(TEST) -d $$f | $(MKDIR) $$f; done)
create_object_directory := $(shell for f in $(OBJECT_DIR); do $(TEST) -d $$f | $(MKDIR) $$f; done)
create_object_library_directory := $(shell for f in $(OBJECT_LIB_DIR); do $(TEST) -d $$f | $(MKDIR) $$f; done)

sources := $(wildcard $(SOURCE_DIR)/*.cpp)
target_objects := $(addprefix $(OBJECT_DIR)/, $(notdir $(TARGET:.cpp=.o)))
lib_objects := $(addprefix $(OBJECT_LIB_DIR)/, $(notdir $(sources:.cpp=.o)))
objects := $(target_objects) $(lib_objects)
dependencies := $(objects:.o=.d)
targets := $(addprefix $(BINARY_DIR)/, $(notdir $(target_objects:.o=)))

# You should indicate whenever a rule does not produce any target output
# with the .PHONY sepcial rule
.PHONY: all
all: $(targets)

# List the prerequisites for building your executable, and fill its
# recipe to tell make what to do with these
$(BINARY_DIR)/%: $(OBJECT_LIB_DIR)/%.o $(lib_objects)
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

# Since your source and object files don't share the same prefix, you
# need to tell make exactly what to do since its built-in rules don't
# cover your specific case
$(OBJECT_LIB_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

#$(OBJECT_DIR)/%.o: $(TARGET_DIR)/%.cpp
#	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

.PHONY: clean
clean:
	$(RM) $(OBJECT_DIR) $(BINARY_DIR)

ifneq "$(MAKECMDGOALS)" "clean"
	-include $(dependencies)
endif

.SECONDARY: $(objects)

.PHONY: help
help:
	@echo 'Build all source files.'
	@echo
	@echo 'Targets in this file are:'
	@echo 'all     Compile and link all source files.'
	@echo 'clean   Remove all intermediate files.'
	@echo 'help    Display this information.'
