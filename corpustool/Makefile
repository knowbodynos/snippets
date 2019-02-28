SHELL := /bin/bash -e

# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export UTIL_DIR	= src/utils
export CXX		= g++
ifeq ($(DEBUG),1)
export CXXFLAGS = -Wall -O0 -g -fno-inline -fkeep-inline-functions -D_FILE_OFFSET_BITS=64 -fPIC -DDEBUG -D_DEBUG
else
export CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC
endif
export LIBS		= -lz \
                  -pthread

PROGRAM_NAME = corpustool

SUBDIRS = $(SRC_DIR)/getcorpus

UTIL_SUBDIRS =	$(UTIL_DIR)/Fasta/ \
				$(UTIL_DIR)/bedFile/ \
				$(UTIL_DIR)/bedFilePE \
				$(UTIL_DIR)/bedFileCorpus \
				$(UTIL_DIR)/gzstream/ \
				$(UTIL_DIR)/fileType/


BUILT_OBJECTS = $(OBJ_DIR)/*.o


all: print_banner $(OBJ_DIR) $(BIN_DIR) $(UTIL_SUBDIRS) $(SUBDIRS)
	@echo "- Building main" $(PROGRAM_NAME) "binary."
	@$(CXX) $(CXXFLAGS) -c src/$(PROGRAM_NAME).cpp -o obj/$(PROGRAM_NAME).o #-I$(UTIL_DIR)/version/
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$(PROGRAM_NAME) $(BUILT_OBJECTS) $(LIBS) #-L$(UTIL_DIR)/BamTools/lib/ -lbamtools
	@echo "done."
.PHONY: all

print_banner:
	@echo "Building " $(PROGRAM_NAME) ":"
	@echo "========================================================="
.PHONY: print_banner

# make the "obj/" and "bin/" directories, if they don't exist
$(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@

# even though these are real directories, treat them as phony targets, forcing to always go in them are re-make.
# a future improvement would be the check for the compiled object, and rebuild only if the source code is newer.
.PHONY: $(UTIL_SUBDIRS) $(SUBDIRS)
$(UTIL_SUBDIRS) $(SUBDIRS): $(OBJ_DIR) $(BIN_DIR)
	@echo "- Building in $@"
	@$(MAKE) --no-print-directory --directory=$@

clean:
	@echo " * Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean