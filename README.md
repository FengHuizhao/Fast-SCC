# Fast-SCC

The source code of paper "Fast-SCC: A Fast Parallel Strongly Connected Components Algorithm via Identifying Fake-link Subgraphs".

## Abstract

Computing strongly connected components (SCC) of a directed graph is among the most fundamental graph-theoretic problems. The most popular parallel SCC detection algorithms are Forward-Backward based on breadth-first search, which are challenged by batch-wise synchronization overhead and redundant operations, especially for the large-diameter graphs on modern multicore platforms. In searching batches, reachability searches may repeatedly inspect adjacency entries that no longer belong to the current unresolved subproblem. In addition, predefined candidate segments may contain resolved vertices or newly exposed 1-SCCs, reducing the amount of useful work performed by a batch.

To tackle these issues, we propose **Fast-SCC**, an optimization of the state-of-the-art (SOTA) GBBS-VGC framework. During local reachability, Fast-SCC identifies fake links from current component labels, narrows per-direction search ranges, and selectively partitions eligible low-degree adjacency lists in place. During pivot selection, it filters resolved and confirmed 1-SCC candidates and uses bounded refill to recover useful frontier positions at controlled cost. These modules reduce redundant work without materializing a separate fake-link graph or changing the SCC membership criterion. We compare Fast-SCC with three SOTA SCC algorithms (GBBS-VGC, GBBS, and Multistep). The results show significant acceleration, with an average speedup of 1.03×, 4.14×, and 2.43×, respectively.

## Repository Structure

```text
Fast-SCC/
├── baseline/
│   ├── GBBS-VGC/
│   │   ├── src/scc               # The executable file of GBBS-VGC
│   │   └── baselines/gbbs/benchmarks/StronglyConnectedComponents/
│   │       └── RandomGreedyBGSS16/StronglyConnectedComponents   # The executable file of GBBS
│   ├── MultiStep/
│   │   └── multistep/scc         # The executable file of MultiStep
├── data/                         # Directed graphs in binary format
├── parlaylib/                    # ParlayLib dependency
├── scripts/
│   ├── ablation.py               # 2 x 2 ablation experiment
│   ├── table4-script.py          # Overall running-time experiment
│   ├── figure8-script.py         # Relative speedup at different thread counts
│   ├── figure9-script.py         # Fast-SCC self-speedup experiment
│   ├── overhead.py               # Module-overhead experiment
│   └── plot/                     # Plotting scripts
└── src/
    ├── scc_paper_swap.cpp
    ├── scc_paper_swap.hpp        # Main Fast-SCC implementation
    ├── scc_paper_swap_ablation.* # Ablation version
    ├── scc_paper_swap_fig10.*    # Search-round instrumentation
    ├── scc_paper_swap_overhead.cpp
    └── Makefile
```

## Requirements

The main implementation requires:

- Linux on an x86-64 processor
- A C++17 compiler, such as GCC 9 or later
- GNU Make
- POSIX threads
- [ParlayLib](https://github.com/cmuparlay/parlaylib)

The experiment scripts additionally require:

- Python 3.9 or later
- Matplotlib
- `numactl` when interleaved memory allocation is enabled
- The baseline implementations used by the corresponding experiment

Install the Python plotting dependency with:

```bash
python3 -m pip install matplotlib
```

## Building

Build all Fast-SCC executables from the repository root:

```bash
make -C src -j
```

The following executables are generated:

| Executable                    | Purpose                                 |
| ----------------------------- | --------------------------------------- |
| `src/scc_paper_swap`          | Main Fast-SCC implementation            |
| `src/scc_paper_swap_ablation` | Four-configuration ablation study       |
| `src/scc_paper_swap_fig10`    | Per-batch search-round instrumentation  |
| `src/scc_paper_swap_overhead` | Fake-link and pivot-selection profiling |

Build only the main implementation with:

```bash
make -C src scc_paper_swap
```

Remove generated binaries and dependency files with:

```bash
make -C src clean
```

The default build uses `-O3`, `-march=native`, and `-mcx16`. To build a portable binary for another compatible x86-64 machine, override the architecture flags:

```bash
make -C src ARCH_FLAGS=
```

A debug build can be generated with:

```bash
make -C src DEBUG=1
```

## Running Fast-SCC

The command format is:

```bash
./src/scc_paper_swap GRAPH [OPTIONS]
```

For example, run Fast-SCC ten times on LiveJournal using 24 worker threads:

```bash
PARLAY_NUM_THREADS=24 \
./src/scc_paper_swap data/soc-LiveJournal1.bin \
  -t 10 -local_reach -local_scc -status
```

The program performs one unmeasured warm-up execution before the measured runs. It reports the running time of each measured execution and their arithmetic mean as `average cost`.

### Main Options

| Option         | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `-t N`         | Execute `N` measured runs; the default is 10                 |
| `-local_reach` | Enable bounded local work in the initial reachability search |
| `-local_scc`   | Enable bounded local work in batched multi-reachability searches |
| `-status`      | Report the number of SCCs and the largest SCC size           |
| `-beta X`      | Set the geometric batch-growth factor; the default is 1.5    |
| `-tau N`       | Override the local-search work budget                        |
| `-lambda N`    | Set the initial hash-bag allocation parameter                |
| `-large`       | Use the large-graph input reader                             |

## Input Graphs

Fast-SCC supports the binary graph representation read by `graph.hpp`. The ten directed graphs used in the experiments:

| Abbreviation | File                    | Download Link                                                |
| ------------ | ----------------------- | ------------------------------------------------------------ |
| LJ           | `soc-LiveJournal1.bin`  | https://www.dropbox.com/scl/fi/836oq1mpruk6y0pul4d4t/soc-LiveJournal1.bin?rlkey=a97bmdoi31a0v87x2j6lmrldg&dl=0 |
| HH5          | `Household.lines_5.bin` | https://www.dropbox.com/scl/fi/ndpp5pq2jjlhuqhpx91r7/Household.lines_5.bin?rlkey=n8s0s9wqyzqfxdl82mmqui1z7&dl=0 |
| CH5          | `CHEM_5.bin`            | https://www.dropbox.com/scl/fi/f0n4o9mrbcfqhw54oa95i/CHEM_5.bin?rlkey=cji03kmi4e5tzmrdj36jsk062&dl=0 |
| GL2          | `GeoLifeNoScale_2.bin`  | https://www.dropbox.com/scl/fi/kckxpp3wcmpbkphikqer4/GeoLifeNoScale_2.bin?rlkey=ybmgdkymbq4i06fvsumk7w9ak&dl=0 |
| GL5          | `GeoLifeNoScale_5.bin`  | https://www.dropbox.com/scl/fi/x3rrl7imfok5742c8g22v/GeoLifeNoScale_5.bin?rlkey=pltytc4h08oyb93hgf3aho0lo&dl=0 |
| GL10         | `GeoLifeNoScale_10.bin` | https://www.dropbox.com/scl/fi/5fw3fhmaqekhihevxvj5e/GeoLifeNoScale_10.bin?rlkey=ed5qn87vdc9m0fjgey7bcopzx&dl=0 |
| GL15         | `GeoLifeNoScale_15.bin` | https://www.dropbox.com/scl/fi/roc5tryp828wj3sdlzceu/GeoLifeNoScale_15.bin?rlkey=20tfzbv69xsh841sobh5tny05&dl=0 |
| GL20         | `GeoLifeNoScale_20.bin` | https://www.dropbox.com/scl/fi/2ufpj603of52ygv8ppqg3/GeoLifeNoScale_20.bin?rlkey=qo9csmyg9bryl3kxejwfhiiwr&dl=0 |
| SQR          | `grid_4000_4000.bin`    | https://www.dropbox.com/scl/fi/7jt1isj9oxejjroeailvk/grid_4000_4000.bin?rlkey=boftsrs25u9gngbmg7mfuaao0&dl=0 |
| REC          | `grid_1000_10000.bin`   | https://www.dropbox.com/scl/fi/wxc3cbdg3i3i6kuqh8ydf/grid_1000_10000.bin?rlkey=ouj3ulurgeomxoruc7g9f6tur&dl=0 |

> [!NOTE]
>
> Graph datasets are **not included** in this repository. Please obtain all `.bin`‑format graph files from [Download Link], place them under the `data/` folder at the project root directory, and then proceed with compilation and execution.

## Experimental Platform

The experiments were conducted on a platform equipped with a 12th Gen Intel Core i9-12900K processor, containing 16 physical cores and 24 hardware threads, with 32 GB of main memory. All compared algorithms were compiled and evaluated under the same environment.

Performance can vary with the compiler, processor architecture, thread placement, memory configuration, and graph structure.

## Citation

If you use this implementation, please cite:

```text
Fast-SCC: A Fast Parallel Strongly Connected Components Algorithm via Identifying Fake-link Subgraphs
```

The complete bibliographic information will be added after publication.

## License

This project is licensed under the MIT License. See `LICENSE` file for details.

