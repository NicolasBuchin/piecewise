CXX=clang++
CC=clang
CXXFLAGS=-std=c++17 -Wall -Wextra
LDFLAGS=-Lblock-aligner/c/target/release -lblock_aligner_c -lstdc++ -Wl,-rpath,$(CURDIR)/block-aligner/c/target/release
INCLUDES=-I.

BALIGNER_SRC = baligner.cpp
BALIGNER_OBJ = $(BALIGNER_SRC:.cpp=.o)

.PHONY: all block_aligner main clean

all: main

block_aligner:
	@echo "Building Rust block-aligner C library..."
	cd block-aligner/c && cargo build --release --features simd_avx2 --offline
	@echo "Generating C header for block-aligner..."
	cd block-aligner/c && cbindgen --config cbindgen.toml --crate block-aligner-c --output ../../block_aligner.h --quiet .

$(BALIGNER_OBJ): $(BALIGNER_SRC) baligner.hpp block_aligner.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

main: block_aligner main.cpp $(BALIGNER_OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o main main.cpp $(BALIGNER_OBJ) $(LDFLAGS)

clean:
	rm -f main $(BALIGNER_OBJ)
	cd block-aligner && cargo clean

