UTILITIES_DIR = ../utils
OBJ_DIR = ../../obj
BIN_DIR = ../../bin

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/Fasta/ \
		   -I$(UTILITIES_DIR)/bedFile/ \
           -I$(UTILITIES_DIR)/bedFileCorpus/ \
           -I$(UTILITIES_DIR)/gzstream/ \
           -I$(UTILITIES_DIR)/lineFileUtilities/ \
           -I$(UTILITIES_DIR)/fileType/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= getcorpusMain.cpp getcorpus.cpp getcorpus.h
OBJECTS= getcorpusMain.o getcorpus.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

# # BUILT_OBJECTS = $(OBJ_DIR)/*.o


# all:
# 	@echo "  * compiling" $(NAME).cpp
# 	@$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/$(NAME).o $(NAME).cpp $(INCLUDES)
# 	# @$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$(NAME) $(BUILT_OBJECTS)
# .PHONY: all

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean