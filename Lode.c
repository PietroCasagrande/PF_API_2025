/*
    Algorithms and Data Structures Project 2024/2025
    Developed by: Pietro Casagrande
    Written with: Kate Text Editor
    Compiled with: GCC 13.3.0 on Ubuntu VM on MacOs

    DEVELOPMENT HISTORY
    - Started with helper functions for init and the construction of the base map via a one-dimensional array, where each element was a "hex" struct of 1 byte cost.
    - Subsequently implemented geometric functions (cube, odd-r offset) to compute distances and neighbors, referencing RedBlobGames theory on hexagonal grids.
    - Added change_cost function to modify hexagon weights according to the formula provided in the specification.
    - Designed air route management with an edge pool and free-list (dynamic allocation only when necessary): memory recycling and AIR_MAX_OUT constraint.
    - Added toggle_air_route and integrated air route management into init and change_cost.
    - travel_cost initially implemented with min_heap (heap_sift_up, heap_sift_down, decrease_key instead of heapify) as priority queue and Dijkstra.
    - First efficient tests: 17.78 s and 7.65 MiB
    - Time optimization with timestamp vectors (sp_dist, sp_done, sp_stamp), reduced by about one second: 16.74 s and 9.54 MiB.
    - Introduction of a cache for travel_cost, implemented via a hashmap with open addressing + linear probing, based on a multiplicative hash with Knuth's constant (inverse golden ratio).
      64-bit Knuth constant for multiplicative hashing: 11400714819323198485ULL (value retrieved from "Fibonacci Hash", programmingpraxis) = 0x9E3779B97F4A7C15 = 2^64 * ((√5 − 1)/2).
      Optimization result: ~10.5 s and ~11.5 MiB.
    - To drop below 10 seconds, introduced neighbor precomputation in nb6 via neighbors_build: result, 8.7s and 26.6MiB.
    - Finally, to get back below 26MiB, reduced the cache size by one order of magnitude.
    - FINAL RESULT: 8.756 s and 25.6 MiB. 
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

/* CONSTANTS DEFINITION */
#define COST_DEFAULT ((cost_t)1)
#define COST_BLOCKED ((cost_t)0)
#define NO_EDGE (-1)                            // Sentinel for "no edge".
#define AIR_MAX_OUT 5                           // Maximum number of air routes per hexagon.
#define DIJK_INF ((int32_t)1e9)                 // "Infinity" for Dijkstra: large enough to avoid overflow when summing weights < 100 many times.
#define TC_CACHE_CAP (1u<<16)                   // 65536 slots, reduced to stay within 26MiB.
#define TC_KNUTH_64 11400714819323198485ULL     // 64-bit Knuth constant for multiplicative hashing: 11400714819323198485ULL = 0x9E3779B97F4A7C15 = 2^64 * ((√5 − 1)/2).

/* TYPE DEFINITIONS */
typedef uint8_t cost_t;             

// hex (hexagon cost)
typedef struct {
    cost_t cost;                   
} hex;

// AirEdge (destination index + route cost + next route index)
typedef struct {
    int32_t to;     
    cost_t  cost;    
    int32_t next;
} AirEdge;

// cube (cubic coordinates q, r, s)
typedef struct {
    int q,r,s;
} cube;  

// directions for cubic coordinates
const cube CUBE_DIRS[6] = {
    {+1, 0, -1},
    {+1, -1, 0},
    {0, -1, +1},
    {-1, 0, +1},
    {-1, +1, 0},
    {0, +1, -1}
};

// Heap_node (index + distance)
typedef struct {
    int32_t node;
    int32_t dist;
} Heap_node;

// Min_heap (elements array + positions array + size + capacity)
typedef struct {
    Heap_node  *arr;    // elements' array
    int32_t    *pos;    // pos[node] = position in the heap (-1 if not there)
    size_t      size;   // how many nodes are actually in the heap
    size_t      cap;    // allocated capacity (max = n_cells)
} Min_heap;

// TC_Entry (src,dst pair compressed into a 64-bit key, the travel cost result dist, current generation)
typedef struct {
    uint64_t    key;
    int32_t     dist;
    int32_t     gen;
} TC_Entry;


/* GLOBAL STATE */

/* Grid */
hex        *map     = NULL;             
size_t      n_cols  = 0;                
size_t      n_rows  = 0;

/* Precomputed neighbors: for each cell, up to 6 neighbors as 1D indices (NO_EDGE if absent) */
int32_t *nb6 = NULL;   // length = n_cells * 6

/* Air Routes*/
int32_t    *air_head        = NULL;         // adjacency list head: one per map cell
AirEdge    *air_edges       = NULL;         // global edge pool
size_t      air_used        = 0;            // number of edges actually in use
size_t      air_cap         = 0;            // allocated capacity
int32_t     air_free_head   = NO_EDGE;      // Free-list of available pool slots: index of first free slot, -1 if empty

Min_heap heap;      // stato globale dell'heap (riutilizzabile)

int32_t *sp_dist = NULL;        // minimum distance found so far for each cell. sp_dist[i] is a reusable array across calls for distances.
uint8_t *sp_done = NULL;        // 0/1 "finalized" node: permanently extracted. sp_done[i] marks nodes already permanently extracted (classic Dijkstra).

int32_t *sp_stamp       = NULL; // timestamp array
int32_t current_stamp   = 0;    // current timestamp

TC_Entry   *tc_tab      = NULL;
uint32_t    tc_cap      = 0;
uint32_t    tc_mask     = 0;
uint32_t    tc_shift    = 0;
int32_t     tc_gen      = 1;


/* FUNCTION PROTOTYPES */
// MAP
void        map_free            (void);                                         // Frees the previous map if one already existed, along with everything that belonged to it.
int         in_bounds           (size_t col, size_t row);                       // Checks whether a hexagon indexed by col, row (offset coordinates) lies within the map boundaries.
size_t      idx                 (size_t col, size_t row);                       // Converts (col, row) of a hexagon to its index in the one-dimensional array (map).

// INIT
void        init                (int columns_in, int rows_in);                  // Initializes the hexagon map with columns_in columns and rows_in rows.

// GEOMETRY - helpers required by change_cost
cube        odd_r_offset_to_cube(size_t col, size_t row);                       // Given the offset coordinates of a hexagon, returns its cubic coordinates.
int         cube_dist           (cube a, cube b);                               // Returns the hexagonal distance (i.e. minimum number of hexagons) between two hexagons given their cubic coordinates.
void        cube_to_odd_r_offset(cube c, int *out_col, int *out_row);           // Given the cubic coordinates of a hexagon, returns its offset coordinates.
int         neighbors_odd_r     (size_t col, size_t row, size_t out_cols[6], size_t out_rows[6]);   // Returns the number of valid neighbors of the central hexagon (passed as input via col, row). The two arrays passed as arguments are "output arrays" containing the coordinates of the neighbors.
int         neighbors_build     (void);                                         // Precomputes neighbors once per init call. Optimization to drop below 10.54 seconds and achieve top marks.

// CHANGE_COST
void        change_cost         (int x_in, int y_in, int v_in, int r_in);

// AIR ROUTES - helpers required by toggle_air_route (pool memory management, logical edge operations, batch utilities)
int         air_reset           (size_t n_cells);                               // Destroys the previous state and re-initializes the air route structures for the new map (one of the functions used in init).
void        air_free_all        (void);                                         // Frees all air route-related structures (used in map_free and air_reset).
int         air_grow            (size_t new_cap);                               // Grows the edge pool (air_edges) when needed, e.g. when a new edge is added.
int32_t     air_alloc_slot      (void);                                         // Allocates a slot for a new edge. First tries to get one from the free-list, otherwise appends growing the pool if full.
void        air_free_slot       (int32_t e);                                    // Recycles a removed slot back into the free-list.

int32_t     air_find_edge       (int32_t src, int32_t dst, int32_t *out_prev);  // Searches the source's adjacency list for the edge "src->dst". Returns the edge index in the pool if found.
int         air_out_degree      (int32_t src);                                  // Counts the number of outgoing air routes from a hexagon (cell). Required to enforce the specification constraint (AIR_MAX_OUT <= 5).
int         air_add_edge        (int32_t src_idx, int32_t dst_idx, cost_t cost);// Adds a new air route to the pool from source to destination, linking it as the head of the source's adjacency list. Uses the free-list or grows the pool.
int         air_remove_edge     (int32_t src_idx, int32_t dst_idx);             // Removes an existing edge from source to destination and recycles the slot into the free-list. Handles both head and middle cases.
cost_t      air_avg_cost_from_src(size_t src_col, size_t src_row);              // Computes the initial cost of a new air route according to the formula provided in the specification.

void        air_apply_delta_from(size_t col, size_t row, int delta);            // Applies a delta (increment or decrement) to all outgoing air routes from the hexagon at (col, row). Used in change_cost to update costs.
void        air_free_from_src   (size_t col, size_t row);                       // Removes all outgoing air routes from hexagon (col, row), recycling each slot into the free-list.

// TOGGLE_AIR_ROUTE
void        toggle_air_route    (int sx_in, int sy_in, int dx_in, int dy_in);

// Heap management functions
int         heap_init           (size_t n_cells);                               // Initializes a heap with maximum capacity equal to the number of map cells.
void        heap_free           (void);                                         // Frees heap-related structures.
void        heap_reset          (void);                                         // Resets the heap without freeing it (globally allocated) at each travel_cost call.
void        heap_swap           (size_t i, size_t j);                           // Swaps nodes in the heap (required by sift_down and sift_up).
void        heap_sift_up        (size_t i);                                     // Moves an element up the heap after a push (insert at the end then bubble up), comparing nodes and maintaining min-heap property (also after decrease_key, which may need to bubble up). Handles the min-heap after insertion.
void        heap_sift_down      (size_t i);                                     // Moves an element down the heap after a pop (place last element at root then sink down comparing with children). Handles the min-heap after minimum removal.
void        heap_push           (int32_t node, int32_t dist);                   // Inserts a new node into the heap.
Heap_node   heap_pop            (void);                                         // Extracts the root of the min-heap (node with minimum distance).
void        heap_decrease_key   (int32_t node, int32_t new_dist);               // Updates the distance of a node already present in the heap.
int         heap_empty          (void);                                         // Checks whether the heap is empty.

// TRAVEL_COST
void        travel_cost         (int sx_in, int sy_in, int dx_in, int dy_in);

// CACHE management functions (optimization for travel_cost)
uint64_t    tc_make_key         (int32_t s, int32_t d);                         // Combines source (shifted left by 32 positions) and destination into a 64-bit key.
uint32_t    tc_index            (uint64_t key);                                 // Computes the hash table index using the multiplication method: m = TC_CACHE_CAP, k = tc_make_key(src, dst), A = TC_KNUTH_64. Floor operation implemented via >>tc_shift.
uint32_t    log2_pow2           (uint32_t x);                                   // Computes the base-two logarithm log2(x) with a loop.
int         tc_init             (void);                                         // Initializes the cache: allocates memory with calloc and computes mask, shift, and gen parameters.
void        tc_free             (void);                                         // Frees cache memory and resets parameters.
void        tc_invalidate       (void);                                         // Invalidates the cache by incrementing the generation counter.
int         tc_get              (int32_t s, int32_t d, int32_t *out_dist);      // Retrieves the distance from the cache (if present), otherwise fails (used in travel_cost, optimizes from O(nlogn) to O(1)).
void        tc_put              (int32_t s, int32_t d, int32_t dist);           // Inserts or updates a (key, distance) pair in the cache.

// MAIN
int main(void) {
    char cmd[32];

    while (scanf("%31s", cmd) == 1) {
        if (strcmp(cmd, "init") == 0) {
            int cols, rows;
            if(scanf("%d %d", &cols, &rows) != 2){
                continue;
            }
            init(cols, rows);
        }
        else if (strcmp(cmd, "change_cost") == 0) {
            int x, y, v, r;
            if(scanf("%d %d %d %d", &x, &y, &v, &r) != 4){
                continue;
            }
            change_cost(x, y, v, r);
        }
        else if (strcmp(cmd, "toggle_air_route") == 0) {
            int col1, row1, col2, row2;
            if (scanf("%d %d %d %d", &col1, &row1, &col2, &row2) != 4){
                continue;
            }
            toggle_air_route(col1, row1, col2, row2);
        }
        else if (strcmp(cmd, "travel_cost") == 0) {
            int col1, row1, col2, row2;
            if(scanf("%d %d %d %d", &col1, &row1, &col2, &row2) != 4){
                continue;
            }
            travel_cost(col1, row1, col2, row2);
        }
    }
    map_free();
    return 0;
}


/* FUNCTION DEFINITIONS in DECLARATION ORDER */

void map_free(void){
    free(map);
    map = NULL;
    n_cols = 0, n_rows = 0;

    air_free_all();
    heap_free();
    free(sp_dist); sp_dist = NULL;
    free(sp_done); sp_done = NULL;

    free(sp_stamp); sp_stamp = NULL;

    free(nb6); nb6 = NULL;  // free precomputed neighbors

    tc_free();              // free the cache
}

int in_bounds(size_t col, size_t row){
    return map != NULL && col < n_cols && row < n_rows;
}

size_t idx(size_t col, size_t row){
    return row * n_cols + col;
}

void init (int columns_in, int rows_in){
    if(columns_in <= 0 || rows_in <= 0) return;     // invalid input

    size_t columns = (size_t) columns_in;
    size_t rows = (size_t) rows_in;

    // cols*rows must fit in size_t 
    if (columns != 0 && rows > SIZE_MAX / columns) return;          // keep the previous map

    size_t n_cells = columns*rows;                                  // n_cells: number of hexagons (or map cells)

    // n_cells*sizeof(hex) must fit in size_t
    if (n_cells != 0 && n_cells > SIZE_MAX / sizeof(hex)) return;   // keep the previous map

    /*freeing the previous map */
    map_free();

    map = (hex*)malloc(n_cells * sizeof(hex));                      // allocate map in memory (one-dimensional array)
    if(!map) return;

    for(size_t i = 0; i < n_cells; ++i) {
        map[i].cost = COST_DEFAULT;                                 // each hexagon's cost initialized to 1
    }

    n_cols = columns;                                               // update to current dimensions
    n_rows = rows;

    if(!air_reset(n_cells)){
        map_free();
        return;
    }

    if(!heap_init(n_cells)) {
        map_free();
        return;
    }

    if(!tc_init()){
        map_free();
        return;
    }

    if (!neighbors_build()){
        map_free();
        return;
    }

    sp_dist = malloc(n_cells * sizeof(int32_t));
    sp_done = malloc(n_cells * sizeof(uint8_t));
    sp_stamp = malloc(n_cells *sizeof(int32_t));
    if(!sp_dist || !sp_done || !sp_stamp) {
        map_free();
        return;
    }

    for(size_t i = 0; i < n_cells; i++) {
        sp_stamp[i] = 0;    // initialize to 0
    }
    current_stamp = 0;

    puts("OK");
}

cube odd_r_offset_to_cube(size_t col, size_t row){
    // conversion from axial coordinates. Source: RedBlobGames pseudocode.
    int r_ax = (int)row;                                        
    int q_ax = (int)col - ((r_ax - (r_ax & 1)) / 2);            
    cube c = (cube){q_ax, r_ax, -q_ax-r_ax};                    
    return c;
}

int cube_dist (cube a, cube b){
    int dq = abs(a.q - b.q);
    int dr = abs(a.r - b.r);
    int ds = abs(a.s - b.s);
    return (dq + dr + ds) / 2;
}

void cube_to_odd_r_offset (cube c, int *out_col, int *out_row){
    int r_ax = c.r;
    int q_ax = c.q;
    int row = r_ax;
    int col = q_ax + ((r_ax - (r_ax & 1)) / 2);     // inverse function of odd_r -> cubic
    *out_col = col;
    *out_row = row;
}

int neighbors_odd_r (size_t col, size_t row, size_t out_cols[6], size_t out_rows[6]) {
    if(!(in_bounds(col, row))) return 0;

    // Cubic coordinates of the central hexagon whose neighbors we want to compute.
    cube C = odd_r_offset_to_cube(col, row);

    int count = 0;
    for(int d = 0; d < 6; ++d) {
        // neighbors in cubic coordinates
        cube N = (cube){
            C.q + CUBE_DIRS[d].q,
            C.r + CUBE_DIRS[d].r,
            C.s + CUBE_DIRS[d].s
        };

        // convert neighbors to offset coordinates
        int n_col_i, n_row_i;
        cube_to_odd_r_offset(N, &n_col_i, &n_row_i);

        // before casting to size_t, check that these neighbors are within map bounds with in_bounds
        if(n_col_i >= 0 && n_row_i >= 0 && (size_t)n_col_i < n_cols && (size_t)n_row_i < n_rows) {
            out_cols[count] = (size_t)n_col_i;
            out_rows[count] = (size_t)n_row_i;
            ++count;
        }
    }
    return count;
}

int neighbors_build(void){
    if (!map || n_cols == 0 || n_rows == 0) return 0;

    size_t n_cells = n_cols * n_rows;
    nb6 = (int32_t*)malloc(n_cells * 6 * sizeof(int32_t));
    if (!nb6) return 0;

    for (size_t r = 0; r < n_rows; ++r){
        for (size_t c = 0; c < n_cols; ++c){
            size_t u_sz = r * n_cols + c;
            int32_t *dst = &nb6[u_sz * 6];

            size_t cols[6], rows[6];
            int k = neighbors_odd_r(c, r, cols, rows);

            /* write valid neighbors */
            for (int i = 0; i < k; ++i){
                dst[i] = (int32_t)(rows[i] * n_cols + cols[i]);
            }
            /* fill remaining slots with NO_EDGE */
            for (int i = k; i < 6; ++i){
                dst[i] = NO_EDGE;
            }
        }
    }
    return 1;
}

void change_cost (int x_in, int y_in, int v_in, int r_in) {
    // Initial validations according to the specification: required to handle invalid input
    if (!map) {puts("KO"); return; }
    if (r_in <= 0) {puts("KO"); return; }
    if (x_in < 0 || y_in < 0) {puts("KO"); return; }
    if (v_in < -10 || v_in > 10) {puts("KO"); return; }

    size_t x = (size_t) x_in;
    size_t y = (size_t) y_in;
    if(!in_bounds(x,y)) {puts("KO"); return; }

    tc_invalidate();

    int v = v_in;   
    int r = r_in;   

    // if v equals 0, the increment is always 0: return the map unchanged
    if (v == 0) {puts("OK"); return; }

    // call function to get the cubic coordinates of the center 
    cube C_center = odd_r_offset_to_cube(x,y);

    // compute the bounding rectangle (square) around the center
    size_t row_min = (y > (size_t)r) ? (y - (size_t)r) : 0;     // if y > (size_t)r, row_min = y - (size_t)r, otherwise row_min = 0
    size_t row_max = y + (size_t)r;
    if (row_max >= n_rows) {
        row_max = n_rows - 1;
    }

    size_t col_min = (x > (size_t)(r+1)) ? (x - (size_t)(r+1)) : 0;     // use (r+1) to account for the column offset on odd rows shifted right
    size_t col_max = x + (size_t)(r+1);
    if (col_max >= n_cols) {
        col_max = n_cols - 1;
    }

    // asymptotic cost is O(r^2)
    for(size_t ry = row_min; ry <= row_max; ++ry) {
        for(size_t cx = col_min; cx <= col_max; ++cx) {

            // Hexagonal distance between the center hexagon and the candidate hexagon, whose cubic coordinates we must compute 
            cube C_cand = odd_r_offset_to_cube(cx, ry);
            int d = cube_dist(C_center, C_cand);

            if (d >= r) continue;   // the specification requires distance to be strictly less than the radius

            /* Increment: floor (v * (r - d) / r ) using INTEGERS ONLY, even for v < 0 */
            int num = r - d;        // >= 0
            int den = r;            // > 0

            int inc;
            if(v >= 0) {
                /* floor( v * num / den) with v, num >= 0: integer division already gives floor */
                int64_t prod = (int64_t)v * (int64_t)num;
                inc = (int)(prod / den);
            } else {
                // case v < 0: floor(v * num / den) = -ceil((-v) * num / den), with ceil(a/b) = (a + b - 1) / b for a >= 0
                // To compute the floor of a negative number, we use: floor(v * num / den) = -ceil((-v) * num / den)
                // To compute ceil(a/b) with positive integers, add (b-1) to the numerator to force rounding up. Example: ceil(5/3) = (5+3-1)/3 = 7/3 = 2 (correct, instead of 5/3=1)
                int64_t a = (int64_t)(-v) * (int64_t)num;   // >= 0
                inc = -(int)((a + den -1) / den);
            }

            // Update + clamp to [0,100]
            int new_cost = (int)map[idx(cx, ry)].cost + inc;
            if(new_cost < 0) new_cost = 0;
            if(new_cost > 100) new_cost = 100;
            map[idx(cx, ry)].cost = (cost_t)new_cost;

            // Apply the same increment "inc" to all outgoing air routes from (cx, ry) 
            if(new_cost == 0) {
                air_free_from_src(cx, ry);
            } else {
                air_apply_delta_from(cx, ry, inc);
            }
        }
    }
    puts("OK");
}

int air_reset(size_t n_cells) {
    // free the previous state
    air_free_all();

    // allocate and set heads to NO_EDGE
    air_head = (int32_t*)malloc(n_cells * sizeof(int32_t));
    if(!air_head) return 0;
    for(size_t i = 0; i < n_cells; ++i) {
        air_head[i] = NO_EDGE;
    }

    // empty pool, will be allocated on demand with air_grow
    return 1;
}

void air_free_all(void){
    free(air_head); air_head = NULL;
    free(air_edges); air_edges = NULL;
    air_used = air_cap = 0;
    air_free_head = NO_EDGE;
}

int air_grow (size_t new_cap) {
    if(new_cap <= air_cap) return 1;

    // Progressive doubling so that it never falls below new_cap 
    size_t target = (air_cap == 0 ? 8 : air_cap*2);
    if(target < new_cap) {
        target = new_cap;
    }

    AirEdge *tmp = (AirEdge*)realloc(air_edges, target * sizeof(AirEdge));
    if (!tmp) return 0;

    air_edges = tmp;
    air_cap = target;
    return 1;
}

int32_t air_alloc_slot(void){
    if(air_free_head != NO_EDGE) {
        // pop from the free-list
        int32_t e = air_free_head;
        air_free_head = air_edges[e].next;
        return e;
    }
    // append at end; if full, grow the pool
    if(air_used == air_cap) {
        if(!air_grow(air_cap == 0 ? 8 : air_cap*2)) {
            return NO_EDGE;
        }
    }
    return (int32_t)air_used++;     // new slot at the end
}

void air_free_slot(int32_t e) {
    air_edges[e].next = air_free_head;
    air_free_head = e;
}

int32_t air_find_edge(int32_t src, int32_t dst, int32_t *out_prev) {
    int32_t prev = NO_EDGE;
    for(int32_t e = air_head[src]; e != NO_EDGE; e = air_edges[e].next) {
        if(air_edges[e].to == dst) {
            if (out_prev){
                *out_prev = prev;
            }
                return e;   // return the found index
        }
        prev = e;
    }
    if (out_prev){
        *out_prev = NO_EDGE;
    }
    return NO_EDGE;
}

int air_out_degree (int32_t src) {
    int deg = 0;
    for(int32_t e = air_head[src]; e != NO_EDGE; e = air_edges[e].next){
        ++deg;
    }
    return deg;
}

int air_add_edge(int32_t src_idx, int32_t dst_idx, cost_t cost) {
    // allocate a slot in the pool
    int32_t e = air_alloc_slot();
    if (e == NO_EDGE) return 0;

    // write the new edge fields (destination and cost)
    air_edges[e].to = dst_idx;
    air_edges[e].cost = cost;

    // link at the head
    air_edges[e].next = air_head[src_idx];
    air_head[src_idx] = e;

    return 1;
}

int air_remove_edge(int32_t src_idx, int32_t dst_idx) {
    int32_t prev = NO_EDGE;

    // Search for edge "e" while tracking the previous one
    int32_t e = air_find_edge(src_idx, dst_idx, &prev);
    if(e == NO_EDGE) return 0;  // edge does not exist

    if(prev == NO_EDGE) {
        // if the removed edge was the list head, move next to head
        air_head[src_idx] = air_edges[e].next;
    } else {
        // if the edge was in the middle, "skip" edge "e"
        air_edges[prev].next = air_edges[e].next;
    }

    // recycle a slot into the free-list
    air_free_slot(e);
    return 1;
}

cost_t air_avg_cost_from_src(size_t src_col, size_t src_row) {
    // one-dimensional index of the source cast to 4 bytes
    size_t src_idx_sz = idx(src_col, src_row);
    int32_t src_idx_i = (int32_t)src_idx_sz;

    // sum and count: include the ground cost of the source hexagon
    int sum = (int)map[src_idx_sz].cost;
    int count = 1;

    // sum of costs of outgoing air routes
    for(int32_t e = air_head[src_idx_i]; e != NO_EDGE; e = air_edges[e].next) {
        sum += (int)air_edges[e].cost;
        ++count;
    }

    int avg = sum / count;      // compute the average
    if(avg < 0) avg = 0;
    if(avg > 100) avg = 100;
    return (cost_t)avg;
}

void air_apply_delta_from(size_t col, size_t row, int delta) {
    if (!map) return;
    if (!in_bounds(col, row)) return;

    int32_t src = (int32_t)idx(col, row);
    for(int32_t e = air_head[src]; e != NO_EDGE; e = air_edges[e].next) {
        int nc = (int)air_edges[e].cost + delta;
        if(nc < 0) nc = 0;
        if(nc > 100) nc = 100;
        air_edges[e].cost = (cost_t)nc;
    }
}

void air_free_from_src(size_t col, size_t row) {
    if (!map) return;
    if (!in_bounds(col, row)) return;

    int32_t src = (int32_t)idx(col, row);
    int32_t e = air_head[src];
    while(e != NO_EDGE) {
        int32_t next = air_edges[e].next;
        air_free_slot(e);                   // free list push
        e = next;
    }
    air_head[src] = NO_EDGE;                // emptied list
}

void toggle_air_route(int sx_in, int sy_in, int dx_in, int dy_in) {
    // Basic validations
    if (!map) {puts("KO"); return; }
    if (sx_in < 0 || sy_in < 0 || dx_in < 0 || dy_in < 0) {puts("KO"); return; }

    size_t sx = (size_t)sx_in;
    size_t sy = (size_t)sy_in;
    size_t dx = (size_t)dx_in;
    size_t dy = (size_t)dy_in;

    if(!in_bounds(sx, sy) || !in_bounds(dx, dy)) {puts("KO"); return; }
    if(!air_head) {puts("KO"); return; }    // air route system not initialized

    // one-dimensional indices for the pool
    size_t src_sz = idx(sx, sy);
    size_t dst_sz = idx(dx, dy);

    // check dimensions fit within int32_t
    if (src_sz > (size_t)INT32_MAX || dst_sz > (size_t)INT32_MAX) {puts("KO"); return; }

    int32_t src = (int32_t)src_sz;
    int32_t dst = (int32_t)dst_sz;

    tc_invalidate();

    // If the route already exists, the command removes it
    if (air_remove_edge(src, dst)) {
        puts("OK");
        return;
    }

    // Otherwise, the command adds it, provided the constraints are satisfied (<= 5)
    if (air_out_degree(src) >= AIR_MAX_OUT) {puts("KO"); return; }

    // Compute the initial cost
    cost_t c_init = air_avg_cost_from_src(sx, sy);

    // Insert (pool and free-list)
    if (!air_add_edge(src, dst, c_init)) {puts("KO"); return; }

    puts("OK");
}

int heap_init (size_t n_cells) {
    heap.arr = (Heap_node*)malloc(n_cells * sizeof(Heap_node));     // allocate maximum space for elements (cap = n_cells), to avoid realloc at runtime
    if(!heap.arr) return 0;
    heap.pos = (int32_t*)malloc(n_cells * sizeof(int32_t));         // allocate space for pos[], the vector tracking current position of node u in the heap (heap.pos[u])
    if(!heap.pos) { free(heap.arr); return 0; }

    heap.size = 0;                                                  
    heap.cap = n_cells;
    for(size_t i = 0; i < n_cells; i++) {
        heap.pos[i] = -1;                                           
    }
    return 1;
}

void heap_free(void) {
    free(heap.arr); heap.arr = NULL;
    free(heap.pos); heap.pos = NULL;
    heap.size = heap.cap = 0;
}

void heap_reset(void) {
    heap.size = 0;
    for(size_t i = 0; i < heap.cap; i++) {
        heap.pos[i] = -1;                                           
    }
}

void heap_swap (size_t i, size_t j) {
    // swap elements at indices i and j 
    Heap_node tmp = heap.arr[i];
    heap.arr[i] = heap.arr[j];
    heap.arr[j] = tmp;

    // update pos[] entries
    heap.pos[heap.arr[i].node] = (int32_t)i;
    heap.pos[heap.arr[j].node] = (int32_t)j;
}

void heap_sift_up (size_t i) {
    while(i > 0) {
        size_t parent = (i - 1) / 2;                            
        if(heap.arr[parent].dist <= heap.arr[i].dist) break;    
        heap_swap(i, parent);                                   
        i = parent;                                             
    }
}

void heap_sift_down (size_t i) {
    while(1) {                                  // sink element at position i until it is smaller than its children
        size_t left = 2*i + 1;                  // left child
        size_t right = 2*i + 2;                 // right child
        size_t smallest = i;

        // find the smallest among i, left and right
        if(left < heap.size && heap.arr[left].dist < heap.arr[smallest].dist) {
            smallest = left;
        }
        if(right < heap.size && heap.arr[right].dist < heap.arr[smallest].dist) {
            smallest = right;
        }

        if(smallest == i) break;                
        heap_swap(i, smallest);                 
        i = smallest;                           
    }
}

void heap_push (int32_t node, int32_t dist) {
    size_t i = heap.size++;                     
    heap.arr[i].node = node;                    
    heap.arr[i].dist = dist;
    heap.pos[node] = (int32_t)i;                
    heap_sift_up(i);                            
}

Heap_node heap_pop (void) {
    Heap_node root = heap.arr[0];               
    heap.pos[root.node] = -1;                   

    heap.size--;                                
    if(heap.size > 0) {                         
        heap.arr[0] = heap.arr[heap.size];      
        heap.pos[heap.arr[0].node] = 0;         
        heap_sift_down(0);                      
    }
    return root;    // return the extracted minimum (node + distance)
}

void heap_decrease_key (int32_t node, int32_t new_dist) {
    int32_t i = heap.pos[node];                
    if (i == -1) return;                       

    heap.arr[i].dist = new_dist;               
    heap_sift_up((size_t)i);                   
}

int heap_empty (void) {
    return heap.size == 0;
}


void travel_cost (int sx_in, int sy_in, int dx_in, int dy_in) {
    if(!map || !air_head || !heap.arr || !sp_dist || !sp_done || !sp_stamp) { puts("-1"); return; }
    if(sx_in < 0 || sy_in < 0 || dx_in < 0 || dy_in < 0) { puts("-1"); return; }

    size_t sx = (size_t)sx_in, sy = (size_t)sy_in, dx = (size_t)dx_in, dy = (size_t)dy_in;
    if(!in_bounds(sx, sy) || !in_bounds(dx, dy)) { puts("-1"); return; }
    if(sx == dx && sy == dy) { puts("0"); return; }

    /*  TIMESTAMP  */
    current_stamp++;                       // new "round"
    // Instead of resetting sp_dist/sp_done at every call, we use timestamps: if sp_stamp[v] != current_stamp, node v has not been touched yet in this Dijkstra round.
    // This avoids O(n) reset operations, replacing them with O(1) checks during the algorithm.
    int32_t s = (int32_t)idx(sx, sy);
    int32_t d = (int32_t)idx(dx, dy);

    int32_t cached;
    if(tc_get(s, d, &cached)) {
        if(cached < 0){
            puts("-1");
        } else {
            printf("%d\n", cached);
        }
        return;
    }

    heap_reset();

    // initialize ONLY the source
    sp_stamp[s] = current_stamp;
    sp_dist[s] = 0;
    sp_done[s] = 0;

    heap_push(s, 0);

    /*  DIJKSTRA  */
    while(!heap_empty()) {
        Heap_node hn = heap_pop();
        int32_t u = hn.node;

        // if u has not been touched in this round (rare but possible), skip
        if(sp_stamp[u] != current_stamp) continue;

        if(sp_done[u]) continue;
        sp_done[u] = 1;
        if(u == d) break;

        cost_t ground = map[u].cost;
        if(ground > 0) {
            // travel by land with precomputed neighbors (final patch for top marks)
            int32_t *vn = &nb6[u * 6];   // pointer to the block of 6 neighbors of u
            for (int i = 0; i < 6; ++i) {
                int32_t v = vn[i];
                if (v == NO_EDGE) continue;

                /* LAZY INIT + relax */
                if (sp_stamp[v] != current_stamp) {
                    sp_stamp[v] = current_stamp;
                    sp_dist[v]  = DIJK_INF;
                    sp_done[v]  = 0;
                }

                int32_t alt = hn.dist + (int32_t)ground;
                if (alt < sp_dist[v]) {
                    sp_dist[v] = alt;
                    if (heap.pos[v] == -1) heap_push(v, alt);
                    else heap_decrease_key(v, alt);
                }
            }

            // travel by air
            for(int32_t e = air_head[u]; e != NO_EDGE; e = air_edges[e].next) {
                int32_t v = air_edges[e].to;

                // **LAZY NEIGHBOR INIT**
                if (sp_stamp[v] != current_stamp) {
                    sp_stamp[v] = current_stamp;
                    sp_dist[v]  = DIJK_INF;
                    sp_done[v]  = 0;
                }

                int32_t alt = hn.dist + (int32_t)air_edges[e].cost;
                if (alt < sp_dist[v]) {
                    sp_dist[v] = alt;
                    if (heap.pos[v] == -1) heap_push(v, alt);
                    else heap_decrease_key(v, alt);
                }
            }
        }
    }

    // result + store in cache
    int32_t result;
    if (sp_stamp[d] != current_stamp || sp_dist[d] == DIJK_INF) result = -1;
    else result = sp_dist[d];

    tc_put(s, d, result);

    if (result < 0) puts("-1");
    else            printf("%d\n", result);
}

uint64_t tc_make_key(int32_t s, int32_t d) {
    return (uint64_t)(uint32_t)s << 32 | (uint64_t)(uint32_t)d; // s in the high bits, d in the low bits
}

// Implementation of the multiplicative method: h(k) = floor(m * frac(k*A)), where A = (√5-1)/2 (golden ratio), m = table size. TC_KNUTH_64 = 2^64 * A, so k*TC_KNUTH_64 gives the fractional part in the high bits.
// Right-shifting by (64-log2(m)) bits is equivalent to multiplying by m and taking floor: this uniformly distributes keys across the hash table
uint32_t tc_index(uint64_t key) {
    uint64_t prod = key * TC_KNUTH_64;      
    return (uint32_t)(prod >> tc_shift);    // floor operation without floats: take the high bits of the product via >>tc_shift (corresponds to floor(m * frac(k*A)) from the course slides)
}

uint32_t log2_pow2(uint32_t x) { // x is a power of 2
    uint32_t p = 0;
    while ((1u << p) != x) ++p;
    return p;
}

int tc_init(void) {
    tc_tab = (TC_Entry*)calloc(TC_CACHE_CAP, sizeof(TC_Entry));
    if(!tc_tab) return 0;
    tc_cap = TC_CACHE_CAP;  // tc_cap is a power of 2 (2^p), so p = log2(tc_cap)
    tc_mask = tc_cap - 1;
    uint32_t p = log2_pow2(tc_cap);     
    // To extract the p most significant bits from a 64-bit product, we shift right by (64-p) positions.
    // For example, if tc_cap = 2^16, p = 16, so to get the top 16 bits as index [0, 2^16-1] we shift by 48 bits.
    tc_shift = 64u - p;
    tc_gen = 1;
    return 1;
}

void tc_free(void){
    free(tc_tab); tc_tab = NULL;
    tc_cap = tc_mask = tc_shift = 0;
    tc_gen = 1;
}

void tc_invalidate(void) {
    if(++tc_gen == INT32_MAX) {
        if(tc_tab) memset(tc_tab, 0, tc_cap * sizeof(TC_Entry));
        tc_gen = 1;
    }
}

int tc_get(int32_t s, int32_t d, int32_t *out_dist){
    if(!tc_tab) return 0;
    uint64_t key = tc_make_key(s,d);
    uint32_t i = tc_index(key);

    while (1){
        TC_Entry *e = &tc_tab[i];
        if (e->gen != tc_gen) return 0;    // "empty slot" for this generation
        if (e->key == key){                // hit
            *out_dist = e->dist;
            return 1;
        }
        i = (i + 1) & tc_mask;             // linear probing
        // To handle collisions, we use linear probing, trying successive slots. Bitwise AND is used instead of division for efficiency.
    }
}

void tc_put(int32_t s, int32_t d, int32_t dist){
    if(!tc_tab) return;
    uint64_t key = tc_make_key(s,d);
    uint32_t i = tc_index(key);

    while (1){
        TC_Entry *e = &tc_tab[i];
        if (e->gen != tc_gen || e->key == key){
            e->key = key;
            e->dist = dist;
            e->gen = tc_gen;
            return;
        }
        i = (i + 1) & tc_mask;
    }
}
