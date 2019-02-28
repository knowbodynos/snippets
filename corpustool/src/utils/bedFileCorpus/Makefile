OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
UTILITIES_DIR = ../../utils/
# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/lineFileUtilities/ \
		   -I$(UTILITIES_DIR)/gzstream/ \
		   -I$(UTILITIES_DIR)/fileType/ \
		   -I$(UTILITIES_DIR)/Fasta/ \
		   -I$(UTILITIES_DIR)/bedFile

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= bedFileCorpus.cpp bedFileCorpus.h
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean