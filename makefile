# Build all source files.
#
# Targets in this file are:
# all     Compile and link all source files.
# clean   Remove all intermediate files.
# help    Display this information.
#
# Copyright (c) 2017 Shota SUGIHARA
# Distributed under the MIT License.

BINARY_DIR := bin
SOURCE_DIR := src
OBJECT_DIR := obj
OBJECT_LIB_DIR := $(OBJECT_DIR)/lib
INCLUDE_DIR := include
TARGET_DIR := src
TARGET := $(wildcard $(TARGET_DIR)/*.cpp)
LIB_DIR := lib

CC := g++
CXX := g++
DEBUGS := -O3
#CXXFLAGS := $(DEBUGS) -std=c++11 -Wall -Wextra -MP -MMD
CXXFLAGS := $(DEBUGS) -std=c++11 -Wextra -MP -MMD -fopenmp -lm
CPPFLAGS := -I $(INCLUDE_DIR) -Iinclude/eigen
LDFLAGS := -L $(LIB_DIR)
#LDLIBS := -lboost_regex
LDLIBS :=

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

.PHONY: all
all: $(targets)

$(BINARY_DIR)/%: $(OBJECT_LIB_DIR)/%.o $(lib_objects)
	$(LINK.cpp) $^ $(LOADLIBES) $(LDLIBS) -o $@

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
	@echo
	@echo 'Copyright (c) 2017 Shota SUGIHARA'
	@echo 'Distributed under the MIT License.'
