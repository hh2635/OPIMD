# Makefile — Exact sampling of Wigner function from open path integral using molecular dynamics simulation.

CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -I../include -I../external

INPUT_DIR    := ../inputs
SAMPLING_DIR := ../sampling
AUTO_DIR     := ../auto
EXEC := main.exe

# 所有源文件
SRC := \
	$(INPUT_DIR)/input.cpp \
	$(SAMPLING_DIR)/OP_init.cpp \
	$(SAMPLING_DIR)/OP_staging.cpp \
	$(SAMPLING_DIR)/OP_force.cpp \
	$(SAMPLING_DIR)/OP_LG.cpp \
	$(AUTO_DIR)/main.cpp

# 让 .o 文件生成在 AUTO_DIR 目录下（即这个 Makefile 所在的目录）
OBJ := $(addprefix $(AUTO_DIR)/, $(notdir $(SRC:.cpp=.o)))

.PHONY: all clean run rebuild

# 链接阶段：使用所有已编译的 .o 文件
all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

$(AUTO_DIR)/main.o: $(AUTO_DIR)/main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(AUTO_DIR)/%.o: $(INPUT_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(AUTO_DIR)/%.o: $(SAMPLING_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(EXEC)
	./$(EXEC)

# 删除可执行文件和本目录下的 .o 文件
clean:
	rm -f $(AUTO_DIR)/*.o $(EXEC)

rebuild: clean all
