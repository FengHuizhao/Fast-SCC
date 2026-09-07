CXX ?= g++

PARLAY_INCLUDE ?= ../parlaylib/include

CPPFLAGS += -I$(PARLAY_INCLUDE)
BASE_FLAGS := -std=c++17 -Wall -Wextra -pthread
ARCH_FLAGS ?= -march=native -mcx16

ifeq ($(DEBUG),1)
OPT_FLAGS := -O0 -g3
else
OPT_FLAGS := -O3
endif

CXXFLAGS ?= $(BASE_FLAGS) $(OPT_FLAGS) $(ARCH_FLAGS)
LDLIBS += -pthread

TARGETS := \
	scc_paper_swap \
	scc_paper_swap_ablation \
	scc_paper_swap_fig10 \
	scc_paper_swap_overhead

COMMON_HEADERS := \
	bfs.hpp \
	get_time.hpp \
	graph.hpp \
	hash_bag.h \
	parseCommandLine.hpp \
	reach.hpp \
	resizable_table.h \
	utilities.h

.PHONY: all main experiments check-parlay clean help

all: $(TARGETS)

main: scc_paper_swap

experiments: \
	scc_paper_swap_ablation \
	scc_paper_swap_fig10 \
	scc_paper_swap_overhead

check-parlay:
	@test -f "$(PARLAY_INCLUDE)/parlay/random.h" || { \
		echo "Error: ParlayLib was not found."; \
		echo "Expected: $(PARLAY_INCLUDE)/parlay/random.h"; \
		echo "Set PARLAY_INCLUDE=/path/to/parlaylib/include if needed."; \
		exit 1; \
	}

scc_paper_swap: scc_paper_swap.cpp scc_paper_swap.hpp $(COMMON_HEADERS) | check-parlay
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

scc_paper_swap_ablation: \
	scc_paper_swap_ablation.cpp \
	scc_paper_swap_ablation.hpp \
	$(COMMON_HEADERS) | check-parlay
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

scc_paper_swap_fig10: \
	scc_paper_swap_fig10.cpp \
	scc_paper_swap_fig10.hpp \
	$(COMMON_HEADERS) | check-parlay
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

scc_paper_swap_overhead: \
	scc_paper_swap_overhead.cpp \
	$(COMMON_HEADERS) | check-parlay
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

clean:
	$(RM) $(TARGETS)

help:
	@echo "Fast-SCC build targets:"
	@echo "  all          Build the main and all instrumented executables"
	@echo "  main         Build only scc_paper_swap"
	@echo "  experiments  Build the ablation, Figure 10, and overhead versions"
	@echo "  clean        Remove generated executables"
	@echo ""
	@echo "Optional variables:"
	@echo "  CXX=<compiler>"
	@echo "  PARLAY_INCLUDE=/path/to/parlaylib/include"
	@echo "  DEBUG=1"
	@echo "  ARCH_FLAGS=   Disable machine-specific compiler flags"