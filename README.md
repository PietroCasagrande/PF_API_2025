# MovHex — Algorithms & Data Structures Final Project
**Politecnico di Milano — A.Y. 2024/2025**
**Grade: 30 cum Laude / 30**

---

## Table of Contents
1. [Project Overview](#project-overview)
2. [Commands & API](#commands--api)
3. [Learning Objectives](#learning-objectives)
4. [Data Structures & Algorithms](#data-structures--algorithms)
5. [Implementation Highlights](#implementation-highlights)
6. [Performance Results](#performance-results)
7. [Development Tools](#development-tools)
8. [Repository Structure](#repository-structure)
9. [Build & Run](#build--run)

---

## Project Overview

**MovHex** is a transport logistics simulation problem. The program models a geographic area as a rectangular grid of **hexagonal tiles**, each with an associated traversal cost. The goal is to efficiently compute minimum-cost routes between hexagons, supporting dynamic updates to costs and the addition/removal of directed air routes between any two hexagons.

### The Map

- The map is a rectangular tiling of hexagons with a fixed number of rows and columns.
- Each hexagon is uniquely identified by its `(column, row)` coordinates (0-indexed, left-to-right, bottom-to-top).
- Every interior hexagon has exactly 6 land neighbors; border hexagons have fewer.
- Each hexagon has an integer traversal cost in the range `[0, 100]`:
  - Cost `0` means the hexagon **cannot be exited** (impassable), but can be a destination.
  - Default cost after `init` is `1`.
- Moving from hexagon A to an adjacent hexagon B costs the traversal cost of A (the **source** hexagon).

### Air Routes

- During execution, directed air routes can be dynamically added or removed between any two hexagons.
- Each hexagon can have **at most 5 outgoing** air routes.
- Using an air route to exit a hexagon costs the **air route's cost** (not the ground cost of the source hexagon).
- The cost of a new air route is the **floor of the average** of all existing outgoing air route costs from the source hexagon and its ground traversal cost:

  ```
  cost = floor( (sum of existing air route costs + ground cost of source) / (num air routes + 1) )
  ```

---

## Commands & API

The program reads commands from `stdin` and writes responses to `stdout`. Each response is terminated by `\n`.

### `init <n_cols> <n_rows>`
Initializes (or re-initializes) a map of `n_rows × n_cols` hexagons. All hexagons get cost `1`. All air routes are cleared.
- **Response:** `OK`

---

### `change_cost <x> <y> <v> <radius>`
Modifies the cost of all hexagons within a circular area (hexagonal distance) centered at `(x, y)` with the given `radius`. The increment applied to each hexagon at hexagonal distance `d` from the center is:

```
increment = floor( v × max(0, (radius - d) / radius) )
```

Costs are clamped to `[0, 100]`. The same delta is also applied to all **outgoing air routes** of each modified hexagon. If a hexagon's cost drops to `0`, all its outgoing air routes are removed.

- **Response:** `OK` on success, `KO` if `(x, y)` is invalid or `radius ≤ 0`.

---

### `toggle_air_route <x1> <y1> <x2> <y2>`
Adds a directed air route from `(x1, y1)` to `(x2, y2)` if it does not exist, or removes it if it does.
- When adding: fails if `(x1, y1)` already has 5 outgoing air routes.
- **Response:** `OK` on success, `KO` otherwise.

---

### `travel_cost <xp> <yp> <xd> <yd>`
Computes the **minimum-cost path** from `(xp, yp)` to `(xd, yd)`.
- **Response:** the minimum cost (integer), `0` if source equals destination, `-1` if unreachable or coordinates invalid.

---

## Learning Objectives

This project is the final exam of the *Algorithms and Principles of Computer Science* course (DEIB, Politecnico di Milano). The goals are:

- Practical application of algorithms and data structures from the course.
- Implementation of an efficient solution to a non-trivial graph problem, paying close attention to **concrete code efficiency** (time and memory).
- Working within strict resource constraints evaluated by an automated grading platform: 6 test batteries with scores `{18, 21, 24, 27, 30, 30L}`. Achieving **30L** requires passing the most demanding battery.

**Constraints enforced by the grader:**
- Language: C (C11, VLA allowed)
- No external libraries beyond the C standard library
- No multithreading
- Input via `stdin`, output via `stdout`
- Compilation flags: `-Wall -Werror -std=gnu11 -O2 -lm`

---

## Data Structures & Algorithms

### Hexagonal Grid — Offset Coordinates & Cube Coordinates

The map uses **odd-r offset coordinates** for indexing, stored as a flat 1D array (`hex *map`), where element `(col, row)` maps to index `row * n_cols + col`. This gives O(1) access by coordinate.

For geometric operations (neighbor computation, hexagonal distance), coordinates are converted to **cube coordinates** `(q, r, s)` with `q + r + s = 0`. The hexagonal distance between two hexagons is:

```c
dist = (|Δq| + |Δr| + |Δs|) / 2
```

This is the standard technique from [RedBlobGames' hexagonal grid theory](https://www.redblobgames.com/grids/hexagons/).

---

### Precomputed Neighbor Table — `nb6`

For each cell, its up to 6 land neighbors are precomputed once at `init` time and stored in the array `nb6` (length = `n_cells × 6`, with sentinel value `NO_EDGE = -1` for missing neighbors). This avoids recomputing geometric neighbor lookups on every Dijkstra iteration, and was the key optimization that pushed runtime below 10 seconds.

---

### Air Route Storage — Edge Pool + Free-List

Air routes are stored using an **adjacency list** model built on a global edge **pool** (`AirEdge *air_edges`):

- `air_head[i]` is the index of the first outgoing air route from cell `i` (or `NO_EDGE`).
- Each `AirEdge` stores `{to, cost, next}` — a singly-linked list per source cell.
- When an edge is removed, its pool slot is **recycled** via a **free-list** (`air_free_head`), so no memory fragmentation occurs and allocation is O(1).
- The pool grows dynamically with progressive doubling (`air_grow`) only when needed.

This design avoids per-edge `malloc`/`free` calls and keeps memory usage predictable.

---

### Dijkstra with Binary Min-Heap

The `travel_cost` function runs **Dijkstra's algorithm** on a mixed graph combining land neighbors and air routes. The priority queue is a custom **binary min-heap** with:

- `heap.arr[]`: array of `(node, dist)` pairs
- `heap.pos[]`: inverse index — `pos[node]` = position of `node` in the heap array (or `-1` if absent)

This `pos` array enables **O(log n) decrease-key** operations (`heap_decrease_key`), which is essential for efficient Dijkstra. The heap is allocated once at `init` and **reset** (not freed) between calls, avoiding repeated memory allocation.

---

### Lazy Initialization with Timestamps — `sp_stamp`

Instead of clearing the distance array `sp_dist[]` and visited array `sp_done[]` at every `travel_cost` call (which would cost O(n)), a **timestamp** mechanism is used:

- `sp_stamp[v]` records the last Dijkstra round in which node `v` was touched.
- `current_stamp` is incremented at each `travel_cost` call.
- If `sp_stamp[v] != current_stamp`, node `v` is treated as uninitialized for this round.

This reduces per-call overhead from O(n) to O(1) per touched node.

---

### Travel Cost Cache — Hash Map with Multiplicative Hashing

Repeated `travel_cost` queries (especially with the same source/destination) are answered in O(1) via a **cache** implemented as a fixed-size open-addressing hash table:

- Key: `(src, dst)` pair packed into a 64-bit integer (`src` in high 32 bits, `dst` in low 32 bits).
- Hash function: **Knuth's multiplicative method** using the 64-bit golden-ratio constant:
  ```
  TC_KNUTH_64 = 11400714819323198485 = 0x9E3779B97F4A7C15 = 2^64 × ((√5 − 1) / 2)
  ```
  ```c
  index = (key * TC_KNUTH_64) >> (64 - log2(capacity))
  ```
- Collision resolution: **linear probing** with bitwise AND masking (table size is always a power of 2).
- **Generational invalidation**: instead of clearing the table on every `change_cost` or `toggle_air_route`, a generation counter `tc_gen` is incremented. Cache entries with a stale generation are treated as empty, making invalidation O(1).

Cache capacity: `2^16 = 65 536` slots (tuned to stay within the 26 MiB memory limit).

---

## Implementation Highlights

| Feature | Technique | Complexity |
|---|---|---|
| Shortest path | Dijkstra + binary min-heap with decrease-key | O((V + E) log V) |
| Neighbor lookup | Precomputed `nb6` table | O(1) per neighbor |
| Dist/visited reset | Timestamp array (`sp_stamp`) | O(1) per node touched |
| Repeated queries | Multiplicative hash map cache | O(1) amortized |
| Air route storage | Pool + free-list singly-linked adjacency lists | O(1) alloc/free |
| Cache invalidation | Generation counter | O(1) |
| `change_cost` sweep | Bounding rectangle + hexagonal distance filter | O(r²) |
| Integer floor division | Avoids floating point for `v < 0` via ceiling identity | O(1) |

---

## Performance Results

The final solution was evaluated on the official DEIB grader (`dum-e.deib.polimi.it`).

**Final score: 8.756 s — 25.6 MiB → 30 cum Laude**

### Optimization History

| Step | Change | Time | Memory |
|---|---|---|---|
| 1 | Baseline Dijkstra + min-heap | 17.78 s | 7.65 MiB |
| 2 | Timestamp vectors (lazy reset) | 16.74 s | 9.54 MiB |
| 3 | Travel cost cache (Knuth multiplicative hash map) | ~10.5 s | ~11.5 MiB |
| 4 | Precomputed neighbor table (`nb6`) | 8.7 s | 26.6 MiB |
| 5 | Reduced cache size (`2^16` slots) | **8.756 s** | **25.6 MiB** |

---

## Development Tools

All development was done on **Ubuntu (VM on macOS)**, using the following toolchain:

### Compilation
```bash
gcc -Wall -Werror -std=gnu11 -O2 -lm Lode.c -o Lode
```
Add `-g3` for debug symbols. The flags mirror those used by the official grader.

### Running Tests
```bash
./Lode < test_input.txt > my_output.txt
diff my_output.txt expected_output.txt
```

### Debugging
| Tool | Purpose | Command |
|---|---|---|
| **GDB** | Runtime state inspection, breakpoints | `gdb ./Lode` |
| **ASAN** | Out-of-bounds, use-after-free detection | compile with `-fsanitize=address` |
| **Memcheck** (Valgrind) | Memory leaks, double-free, uninitialized reads | `valgrind ./Lode < input.txt` |

### Profiling
| Tool | Purpose | Command |
|---|---|---|
| **Callgrind** (Valgrind) | CPU time profiling, call graph | `valgrind --tool=callgrind ./Lode < input.txt` then `kcachegrind callgrind.out.PID` |
| **Massif** (Valgrind) | Heap memory usage over time | `valgrind --tool=massif ./Lode < input.txt` then `massif-visualizer massif.out.PID` |
| **`/usr/bin/time`** | Peak memory (RSS) summary | `/usr/bin/time -v ./Lode < input.txt` |

### Install all tools (Debian/Ubuntu)
```bash
sudo apt install gdb valgrind kcachegrind massif-visualizer build-essential
```

---

## Repository Structure

```
.
├── Lode.c                         # Final submission — 30L solution
├── README.md                      # This file
│
├── doc/                           # Reference documentation (Italian)
│   ├── specifica_2024_2025.pdf    # Official project specification
│   ├── Slides_PFAPI_24_25.pdf     # Lecture slides (problem description)
│   └── strumenti_progetto_api.pdf # Tools guide (GDB, Valgrind, ASAN, compilation)
│
├── testing/
│   ├── generators/                # Automated test case generators
│   │   ├── generator_mac          # Input generator binary (macOS)
│   │   └── solver_mac             # Reference solver binary (macOS)
│   └── public_tests/              # Public test cases provided by the course
│       ├── edge_cases.txt
│       ├── edge_cases.result
│       └── ...
│
└── previous_submissions/          # Earlier submissions with lower scores
    └── 24.c                       # Submission that achieved 24/30
```

---

## Build & Run

```bash
# Compile
gcc -Wall -Werror -std=gnu11 -O2 -lm Lode.c -o Lode

# Run with a test file
./Lode < testing/public_tests/edge_cases.txt

# Run and compare output
./Lode < testing/public_tests/edge_cases.txt > out.txt
diff out.txt testing/public_tests/edge_cases.result
```

---

*Project developed by Pietro Casagrande — DEIB, Politecnico di Milano, A.Y. 2024/2025.*
