#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

/* DEFINIZIONE COSTANTI */
#define COST_DEFAULT ((cost_t)1)
#define COST_BLOCKED ((cost_t)0)
#define NO_EDGE (-1)    // Sentinella per "nessun arco"
#define AIR_MAX_OUT 5   // massimo numero di rotte aeree per esagono
#define INF32 0x3f3f3f3f        // no overflow quando sommo pesi < 100 molte volte


/* DEFINIZIONE TIPI */
typedef uint8_t cost_t;             // 0...100

typedef struct {
    cost_t cost;                    // the only member of the struct hex is the cost, which occupies 2 bytes
} hex;

typedef struct {
    int32_t to;     // Indice 1d della destinazione (0, ..., n_cells-1), oppure NO_EDGE
    cost_t cost;    // costo della rotta (clamp 0...100)
    int32_t next;   // indice del prossimo arco nella lista (stessa sorgente), oppure NO_EDGE
} AirEdge;

typedef struct {int q,r,s;} cube;   /* cubic coordinates (q,r,s) with the constraint q+r+s = 0 */

const cube CUBE_DIRS[6] = {
    {+1, 0, -1},
    {+1, -1, 0},
    {0, -1, +1},
    {-1, 0, +1},
    {-1, +1, 0},
    {0, +1, -1}
};

// Heap node
typedef struct {
    int32_t node;   // cell's index
    int32_t dist;   // current distance from the source
} Heap_node;

// Heap structure
typedef struct {
    Heap_node  *arr;    // elements' array
    int32_t    *pos;    // pos[node] = position in the heap (-1 if not there)
    size_t      size;   // how many nodes are actually in the heap
    size_t      cap;    // allocated capacity (max = n_cells)
} Min_heap;


/* STATO GLOBALE */

/* Grid */
hex        *map     = NULL;             // map is a 1d array where each cell is a struct hex
size_t      n_cols  = 0;                // the number of columns and rows are initialized at zero
size_t      n_rows  = 0;

/* Air Routes*/
int32_t    *air_head        = NULL;         // testa della lista: una per ogni cella della mappa
AirEdge    *air_edges       = NULL;         // pool globale di archi
size_t      air_used        = 0;            // quanti archi realmente usati
size_t      air_cap         = 0;            // capacità allocata
int32_t     air_free_head   = NO_EDGE;      // Free-list degli slot liberi nel pool: indice del primo libero, -1 se vuota


Min_heap heap;      // heap's global state (reusable)

int32_t *sp_dist = NULL;        // distanza minima trovata finora per ogni cella. L'array sp_dist[i] è un array riusabile tra chiamate per le distanze.
uint8_t *sp_done = NULL;        // 0/1 nodo "finalizzato": estratto in modo definitivo. L'array sp_done[i] marca i nodi già estratti in modo definitivo (classico Dijkstra).


/* PROTOTIPI */
// MAPPA
void map_free(void);
static inline int in_bounds(size_t col, size_t row);
static inline size_t idx(size_t col, size_t row);

// INIT
void init (int columns_in, int rows_in);

// GEOMETRIA - helper necessari a change_cost
static inline cube odd_r_offset_to_cube(size_t col, size_t row);
static inline int cube_dist (cube a, cube b);
// INUTILIZZATA: static inline int hex_distance_odd_r_offset (size_t col1, size_t row1, size_t col2, size_t row2);
void cube_to_odd_r_offset (cube c, int *out_col, int *out_row);
int neighbors_odd_r (size_t col, size_t row, size_t out_cols[6], size_t out_rows[6]);

// CHANGE_COST
void change_cost (int x_in, int y_in, int v_in, int r_in);

// ROTTE AEREE - helper necessari a toggle_air_route
// gestione memoria pool
int air_reset(size_t n_cells);
void air_free_all(void);
int air_grow (size_t new_cap);
int32_t air_alloc_slot(void);
void air_free_slot(int32_t e);

// operazioni logiche sugli archi
int32_t air_find_edge(int32_t src, int32_t dst, int32_t *out_prev);
int air_out_degree (int32_t src);
int air_add_edge(int32_t src_idx, int32_t dst_idx, cost_t cost);
int air_remove_edge(int32_t src_idx, int32_t dst_idx);
cost_t air_avg_cost_from_src(size_t src_col, size_t src_row);

// batch utili
void air_apply_delta_from(size_t col, size_t row, int delta);
void air_free_from_src(size_t col, size_t row);

// TOGGLE_AIR_ROUTE
void toggle_air_route(int sx_in, int sy_in, int dx_in, int dy_in);

// TODO: TRAVEL_COST

/* SHORTEST PATH wih DIJKSTRA'S ALGORITHM implemented trough PRIORITY QUEUES, realized with MIN-HEAP
 * Used as a priority queue in Dijkstra's algorithm.
 * The key is the shortest distance (int32_t, clamped in [0..100])
 * Each node represents a cell (an hexagon of the grid)
 * The heap supports:
 *      - push (insert)
 *      - pop (minimum estraction)
 *      - decrease_key (update best (shortest) distance)
 * The min-heap implemented is reusable, as it is allocated once inside init(), and later resetted
 */
// Heap functions
int         heap_init           (size_t n_cells);                               // initializes the heap with max capacity equal to the number of cells (heap.cap = n_cells)
void        heap_free           (void);                                         // frees the structures of the heap
void        heap_reset          (void);                                         // heap reset without free, used before each travel_cost
void        heap_swap           (size_t i, size_t j);                           // swap helper
void        heap_sift_up        (size_t i);                                     // "climbs" back the heap until (dist) property is respected
void        heap_sift_down      (size_t i);                                     // "climbs" down the heap
void        heap_push           (int32_t node, int32_t dist);                   // adds a new node in the heap
Heap_node   heap_pop            (void);                                         // extracts the node with minimum distance
void        heap_decrease_key   (int32_t node, int32_t new_dist);               // updates the distance of a node which is already inside the heap
int         heap_empty          (void);                                         // checks if the heap is empty
void        travel_cost         (int sx_in, int sy_in, int dx_in, int dy_in);


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
                continue;   //input malformato
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
        else if (strcmp(cmd, "print") == 0) {
            for (size_t r=0; r<n_rows; r++) {
                if (r & 1) printf("  ");        // indentazione per odd-r
                for (size_t c=0; c<n_cols; c++) printf("%3u ", map[idx(c,r)].cost);
                puts("");
            }
        }
        /*
        else if (strcmp(cmd, "print_air_from") == 0) {
            int cx, cy; if(scanf("%d %d", &cx, &cy) != 2) continue;
            if (cx >= 0 && cy >= 0 && in_bounds((size_t)cx,(size_t)cy) && air_head) {
                int32_t s = (int32_t)idx((size_t)cx,(size_t)cy);
                printf("Air from (%d,%d):", cx, cy);
                for (int32_t e = air_head[s]; e != NO_EDGE; e = air_edges[e].next) {
                    int32_t dst = air_edges[e].to;
                    size_t drow = (size_t)(dst / n_cols);
                    size_t dcol = (size_t)(dst % n_cols);
                    printf(" -> (%zu,%zu)[%u]", dcol, drow, (unsigned)air_edges[e].cost);
                }
                puts("");
            } else {
                puts("KO");
            }
        }

        else if (strcmp(cmd, "exit") == 0) {
            map_free();
	    break;
        }
        else {
            printf("Comando non riconosciuto: %s\n", cmd);
        }
        */
    }
    map_free();
    return 0;
}


/* DEFINIZIONI in ORDINE di DICHIARAZIONE*/

/*frees the previous map in case it already exists*/
void map_free(void){
    free(map);
    map = NULL;
    n_cols = 0, n_rows = 0;

    air_free_all();
    heap_free();
    free(sp_dist); sp_dist = NULL;
    free(sp_done); sp_done = NULL;
}

// IN_BOUNDS: Returns 1 (true) if (col,row) is inside the map and the map exists
static inline int in_bounds(size_t col, size_t row){
    return map != NULL && col < n_cols && row < n_rows;
}

// IDX: Converts (col, row) into a index idx of the 1d array. Useful only if in_bounds is TRUE
static inline size_t idx(size_t col, size_t row){
    return row * n_cols + col;
}

/*INIT: the init function takes the number of columns and the number of rows as parameters, and creates a grid of columns*rows dimension, which is stored in a 1 dimension dynamic array named as "*map", through a malloc call that gives to the pointer the right size, which is the number of cells of the array (columns*rows), times the sizeof the struct hex (which is just 2 bytes of the member cost)*/
void init (int columns_in, int rows_in){
    if(columns_in <= 0 || rows_in <= 0) return;     // invalid input

    size_t columns = (size_t) columns_in;
    size_t rows = (size_t) rows_in;

    /* cols*rows deve stare in size_t */
    if (columns != 0 && rows > SIZE_MAX / columns) return;          // no OK, the previous map remains

    size_t n_cells = columns*rows;                                  // n_cells: number of cells

    /* n_cells*sizeof(hex) deve stare in size_t */
    if (n_cells != 0 && n_cells > SIZE_MAX / sizeof(hex)) return;   // no OK, the previous map remains

    /*freeing the previous map */
    map_free();

    map = (hex*)malloc(n_cells * sizeof(hex));                      // MEMORY ALLOCATION for the 1d ARRAY
    if(!map) return;

    for(size_t i = 0; i < n_cells; ++i) {
        map[i].cost = COST_DEFAULT;                                 // each exagon's cost init to 1
    }

    n_cols = columns;                                               // update of current dimension
    n_rows = rows;

    if(!air_reset(n_cells)){
        map_free();
        return;
    }

    if(!heap_init(n_cells)) {
        map_free();
        return;
    }

    sp_dist = malloc(n_cells * sizeof(int32_t));
    sp_done = malloc(n_cells * sizeof(uint8_t));
    if(!sp_dist || !sp_done) {
        map_free();
        return;
    }

    puts("OK");
}


/* GEOMETRY FUNCTIONS */

/* odd_r_offset_to_cube is the function that given the offset coordinates of an hex, returns its cubic coordinate */
static inline cube odd_r_offset_to_cube(size_t col, size_t row){
    int r_ax = (int)row;                                        // "axial"-r = row
    int q_ax = (int)col - ((r_ax - (r_ax & 1)) / 2);            // "axial"-q is given by the odd-r to cubic formula
    cube c = (cube){q_ax, r_ax, -q_ax-r_ax};                    // (q,r,s) where s = -q-r in axial coordinates
    return c;
}

/* cube_dist: returns the hexagonal distance between two hexagons given their cubic coordinates */
static inline int cube_dist (cube a, cube b){
    int dq = abs(a.q - b.q);
    int dr = abs(a.r - b.r);
    int ds = abs(a.s - b.s);
    return (dq + dr + ds) / 2;
}

/* hex_distance_odd_r_offset: returns the distance between two hexagons given their offset coordinates 
static inline int hex_distance_odd_r_offset (size_t col1, size_t row1, size_t col2, size_t row2){
    cube a = odd_r_offset_to_cube(col1, row1);
    cube b = odd_r_offset_to_cube(col2, row2);
    return cube_dist(a, b);
} */ // E' INUTILIZZATA

/* CUBE_TO_ODD_R_OFFSET: it is the inverse function of ODD_R_OFFSET_TO_CUBE, as it returns the offset coordinates of an hex, given its cubic coordinates. */
// We are passing two pointers in order to be given in return our col and row separatly.
void cube_to_odd_r_offset (cube c, int *out_col, int *out_row){
    int r_ax = c.r;
    int q_ax = c.q;
    int row = r_ax;
    int col = q_ax + ((r_ax - (r_ax & 1)) / 2);     // inverse function of odd_r -> cubic
    *out_col = col;
    *out_row = row;
}

/* NEIGHBORS_ODD_R : La funzione ritorna il numero di vicini validi dell'esagono centrale. I due array passati come argomenti sono "array di output", ovvero contengono le coordinate dei vicini. */
int neighbors_odd_r (size_t col, size_t row, size_t out_cols[6], size_t out_rows[6]) {
    if(!(in_bounds(col, row))) return 0;

    /* Center in cubic coordinates */
    cube C = odd_r_offset_to_cube(col, row);

    int count = 0;
    for(int d = 0; d < 6; ++d) {
        // neighbors in cubic coordinates: adding the direction d
        cube N = (cube){
            C.q + CUBE_DIRS[d].q,
            C.r + CUBE_DIRS[d].r,
            C.s + CUBE_DIRS[d].s
        };

        // conversion of the neighbors coordinates to offset (return type is int to manage negatives)
        int n_col_i, n_row_i;
        cube_to_odd_r_offset(N, &n_col_i, &n_row_i);

        // before the cast to size_t, we'll do a bounds check
        if(n_col_i >= 0 && n_row_i >= 0 && (size_t)n_col_i < n_cols && (size_t)n_row_i < n_rows) {
            out_cols[count] = (size_t)n_col_i;
            out_rows[count] = (size_t)n_row_i;
            ++count;
        }
    }
    return count;
}

/* CHANGE_COST */
void change_cost (int x_in, int y_in, int v_in, int r_in) {
    /* Validazioni iniziali secondo la specifica: necessarie per gestire il caso di input errato */
    if (!map) {puts("KO"); return; }
    if (r_in <= 0) {puts("KO"); return; }
    if (x_in < 0 || y_in < 0) {puts("KO"); return; }
    if (v_in < -10 || v_in > 10) {puts("KO"); return; }

    size_t x = (size_t) x_in;
    size_t y = (size_t) y_in;
    if(!in_bounds(x,y)) {puts("KO"); return; }

    int v = v_in;   // can be negative
    int r = r_in;   // garanteed to be greater than zero

    // if v = 0, the increment is always 0, therefore we respect the spec and we output OK
    if (v == 0) {puts("OK"); return; }

    // getting the cubic coordinates of the center
    cube C_center = odd_r_offset_to_cube(x,y);

    // we'll now calculate a rectangle around the center
    size_t row_min = (y > (size_t)r) ? (y - (size_t)r) : 0;     // if the expression y > (size_t)r is true, row_min = y - (size_t)r, else row_min = 0;
    size_t row_max = y + (size_t)r;
    if (row_max >= n_rows) {
        row_max = n_rows - 1;
    }

    size_t col_min = (x > (size_t)(r+1)) ? (x - (size_t)(r+1)) : 0;     // we're using (r+1) to catch the offset on odd-r columns
    size_t col_max = x + (size_t)(r+1);
    if (col_max >= n_cols) {
        col_max = n_cols - 1;
    }

    // asyntothic cost O(r^2)
    for(size_t ry = row_min; ry <= row_max; ++ry) {
        for(size_t cx = col_min; cx <= col_max; ++cx) {

            /* Distanza esagonale tra l'esagono al centro del rettangolo e l'esagono candidato, di cui dobbiamo calcolare le coordinate cubiche */
            cube C_cand = odd_r_offset_to_cube(cx, ry);
            int d = cube_dist(C_center, C_cand);

            if (d >= r) continue;   // la specifica richiede che la distanza sia strettamente minore del raggio

            /* Incremento: floor (v * (r - d) / r ) con SOLO interi, anche per v < 0 */
            int num = r - d;        // >= 0
            int den = r;            // > 0

            int inc;
            if(v >= 0) {
                /* floor( v * num / den) con v, num >= 0: la / intera fa già il floor */
                int64_t prod = (int64_t)v * (int64_t)num;
                inc = (int)(prod / den);
            } else {
                // caso v < 0: floor(v * num / den) = - ceil ( (-v) * num / den), con ceil(a/b) = (a + b -1) / b per a >= 0
                int64_t a = (int64_t)(-v) * (int64_t)num;   // >= 0
                inc = -(int)((a + den -1) / den);
            }

            /* Aggiornamento + clamp in [0,100] */
            int new_cost = (int)map[idx(cx, ry)].cost + inc;
            if(new_cost < 0) new_cost = 0;
            if(new_cost > 100) new_cost = 100;
            map[idx(cx, ry)].cost = (cost_t)new_cost;

            /* Applicare lo stesso incremento "inc" alle rotte aeree uscenti di (cx, ry) */
            if(new_cost == 0) {
                air_free_from_src(cx, ry);
            } else {
                air_apply_delta_from(cx, ry, inc);
            }
        }
    }
    puts("OK");
}



// POOL'S MEMORY MENAGEMENT

// 11) AIR_RESET: "destroys" the previous state and re initialize the structs for the new map, used in INIT()
int air_reset(size_t n_cells) {
    // "destruction" of the previous state
    air_free_all();

    // alloca e setta teste a NO_EDGE
    air_head = (int32_t*)malloc(n_cells * sizeof(int32_t));
    if(!air_head) return 0;
    for(size_t i = 0; i < n_cells; ++i) {
        air_head[i] = NO_EDGE;
    }

    // pool vuoto, verra allocato on demand con air_grow
    return 1;
}

// 12) AIR_FREE_ALL: frees all the structures of the air routes, used in MAP_FREE() and at the beginning of AIR_RESET() as well
void air_free_all(void){
    free(air_head); air_head = NULL;
    free(air_edges); air_edges = NULL;
    air_used = air_cap = 0;
    air_free_head = NO_EDGE;
}

// 1) AIR_GROW: grows the pool "air_edges[]" when needed. It is used in AIR_ALLOC_SLOT, and therefore, indirectly, in AIR_ADD_EDGE
int air_grow (size_t new_cap) {
    if(new_cap <= air_cap) return 1;

    // progressive doubling, never less than new_cap
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

// 2) AIR_ALLOC_SLOT: obtains a slot for a new edge. It first tries from the free-list, otherwise appends growing the pool if full.
// It returns a free index from the pool if there is one.
int32_t air_alloc_slot(void){
    if(air_free_head != NO_EDGE) {
        // pop from the free-list
        int32_t e = air_free_head;
        air_free_head = air_edges[e].next;
        return e;
    }
    // tail append, growing if full
    if(air_used == air_cap) {
        if(!air_grow(air_cap == 0 ? 8 : air_cap*2)) {
            return NO_EDGE;
        }
    }
    return (int32_t)air_used++;     // new slot in the tail
}

// 3) AIR_FREE_SLOT: recycles a removed slot from the free list
void air_free_slot(int32_t e) {
    air_edges[e].next = air_free_head;
    air_free_head = e;
}


// LOGIC OPERATIONS on the POOL

// 4) AIR_FIND_EDGE: searchs in the adjacency list of src the edge src->dst. Returns the edge's index in the pool if found, else NO_EDGE.
// If out_prev != NULL, the function writes the index of the previous node in the chain (NO_EDGE if it was the head)
int32_t air_find_edge(int32_t src, int32_t dst, int32_t *out_prev) {
    int32_t prev = NO_EDGE;
    for(int32_t e = air_head[src]; e != NO_EDGE; e = air_edges[e].next) {
        if(air_edges[e].to == dst) {
            if (out_prev){
                *out_prev = prev;
            }
                return e;   // ritorna l'indice trovato
        }
        prev = e;
    }
    if (out_prev){
        *out_prev = NO_EDGE;
    }
    return NO_EDGE;
}

// 5) AIR_OUT_DEGREE: counts how many air routes does the src cell have
// It is needed in order to respect the bond of "max 5" air routes
int air_out_degree (int32_t src) {
    int deg = 0;
    for(int32_t e = air_head[src]; e != NO_EDGE; e = air_edges[e].next){
        ++deg;
    }
    return deg;
}

// 7) AIR_ADD_EDGE: adds a new air route from src to dst in the pool, linking it as the head of the src's list
// It uses the free list if there are any free slots, otherwise it grews the pool and appends.
int air_add_edge(int32_t src_idx, int32_t dst_idx, cost_t cost) {
    // allocates a slot in the pool (free-list or append)
    int32_t e = air_alloc_slot();
    if (e == NO_EDGE) return 0;     // failed re-allocation

    // writes the new members of the edge
    air_edges[e].to = dst_idx;
    air_edges[e].cost = cost;

    // head linking: next = current head; new head = e;
    air_edges[e].next = air_head[src_idx];
    air_head[src_idx] = e;

    return 1;
}

// 8) AIR_REMOVE_EDGE: removes an already existing edge from src to dst from the "chain" of src, and recycles the slot in free-list.
// The function handles both the case in which the edge is in the head, both the case in which it is in the middle.
int air_remove_edge(int32_t src_idx, int32_t dst_idx) {
    int32_t prev = NO_EDGE;

    // Search the edge "e", remembering the previous one
    int32_t e = air_find_edge(src_idx, dst_idx, &prev);
    if(e == NO_EDGE) return 0;  // the edge doesn't exist

    if(prev == NO_EDGE) {
        // the removed edge was the head of the list: moves head to next
        air_head[src_idx] = air_edges[e].next;
    } else {
        // removed edge was in the middle: skips "e"
        air_edges[prev].next = air_edges[e].next;
    }

    // recycles the slot in free-list
    air_free_slot(e);
    return 1;
}

/* 6) AIR_AVG_COST_FROM_SRC: calculates the beginning cost of a new air route from src as the mean (floor) between the ground-cost of the src hex and
 * the costs of all the already existing air routes exiting from that hex. If there are no air routes exiting from the hex, the cost will e the ground
 * one.
 */
cost_t air_avg_cost_from_src(size_t src_col, size_t src_row) {
    // 1d index of the src, both in size_t and int32_t versions
    size_t src_idx_sz = idx(src_col, src_row);
    int32_t src_idx_i = (int32_t)src_idx_sz;

    // sum and count: we always include the ground_cost of the src hexagon
    int sum = (int)map[src_idx_sz].cost;
    int count = 1;

    // sums the exiting air routes costs
    for(int32_t e = air_head[src_idx_i]; e != NO_EDGE; e = air_edges[e].next) {
        sum += (int)air_edges[e].cost;
        ++count;
    }

    int avg = sum / count;      // integers' floor
    if(avg < 0) avg = 0;
    if(avg > 100) avg = 100;
    return (cost_t)avg;
}


// USEFUL BATCH OPERATIONS

// 9) AIR_APPLY_DELTA_FROM: it applies an increment (or decrement) delta to all the exiting air routes from (col, row), with [0..100] clamp.
// This function will be used inside change_cost, shortly after the hex's cost update.
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

// 10) AIR_FREE_FROM_SRC: removes evey air route exiting from (col, row), recycling each slot in the free-list.
// Function that will be used later.
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


/* TOGGLE_AIR_ROUTE:
 * Input: (sx, sy) source, (dx, dy) destination (odd-r offset coordinates)
 * If the route (src->dst) already exists, the function removes it and prints "OK"
 * Otherwise, it adds the route:
 *      - it respects the max number of air routes exiting from an hex AIR_MAX_OUT
 *      - initial cost = floor ( mean ( ground_cost(src) + exiting_air_routes_costs(src))), clamped in [0..100]
 *      - head insert: prints "OK"
 * KO is for invalid input | out of bounds | failed alloc | not ready routes system
 */
void toggle_air_route(int sx_in, int sy_in, int dx_in, int dy_in) {
    // Base validations
    if (!map) {puts("KO"); return; }
    if (sx_in < 0 || sy_in < 0 || dx_in < 0 || dy_in < 0) {puts("KO"); return; }

    size_t sx = (size_t)sx_in;
    size_t sy = (size_t)sy_in;
    size_t dx = (size_t)dx_in;
    size_t dy = (size_t)dy_in;

    if(!in_bounds(sx, sy) || !in_bounds(dx, dy)) {puts("KO"); return; }
    if(!air_head) {puts("KO"); return; }    // unitialized routes system

    // 1d indixes for pool
    size_t src_sz = idx(sx, sy);
    size_t dst_sz = idx(dx, dy);

    // check size (inside int32_t)
    if (src_sz > (size_t)INT32_MAX || dst_sz > (size_t)INT32_MAX) {puts("KO"); return; }

    int32_t src = (int32_t)src_sz;
    int32_t dst = (int32_t)dst_sz;

    // If the route already exists, remove it
    if (air_remove_edge(src, dst)) {
        puts("OK");
        return;
    }

    // Else, add the new route, but be sure to respect the max number
    if (air_out_degree(src) >= AIR_MAX_OUT) {puts("KO"); return; }

    // Setting of the initial cost
    cost_t c_init = air_avg_cost_from_src(sx, sy);

    // Insert (pool and free-list)
    if (!air_add_edge(src, dst, c_init)) {puts("KO"); return; }

    puts("OK");
}



// HEAP_INIT: initializes the heap with max capacity equal to the number of cells (heap.cap = n_cells)
int heap_init (size_t n_cells) {
    heap.arr = (Heap_node*)malloc(n_cells * sizeof(Heap_node));     // alloca lo spazio massimo per gli elementi (cap = n_cells), cosi da evitare realloc in run-time
    if(!heap.arr) return 0;
    heap.pos = (int32_t*)malloc(n_cells * sizeof(int32_t));         // alloca lo spazio per pos[], il vettore che tiene traccia della posizione corrente del nodo u in heap (heap.pos[u])
    if(!heap.pos) { free(heap.arr); return 0; }

    heap.size = 0;                                                  // size = 0 -> l'heap è vuoto
    heap.cap = n_cells;
    for(size_t i = 0; i < n_cells; i++) {
        heap.pos[i] = -1;                                           // no nodes inside the heap at the initialization
    }
    return 1;
}

// HEAP_FREE
void heap_free(void) {
    free(heap.arr); heap.arr = NULL;
    free(heap.pos); heap.pos = NULL;
    heap.size = heap.cap = 0;
}

// HEAP_RESET: the heap is resetted without using free. This function is used before each travel_cost
void heap_reset(void) {
    heap.size = 0;
    for(size_t i = 0; i < heap.cap; i++) {
        heap.pos[i] = -1;                                           // no nodes inside the heap
    }
}

// HEAP_SWAP: swap helper needed to swap the position between two nodes
void heap_swap (size_t i, size_t j) {
    // swap between the elements at the indexes i and j
    Heap_node tmp = heap.arr[i];
    heap.arr[i] = heap.arr[j];
    heap.arr[j] = tmp;

    // updates the entries of pos[] for the swapped nodes, so that pos[u] always makes sense
    heap.pos[heap.arr[i].node] = (int32_t)i;
    heap.pos[heap.arr[j].node] = (int32_t)j;
}

// HEAP_SIFT_UP: climbs up the heap until dist is respected (the order of the min heap is respected)
// Si usa dopo un push (inserimento in coda e poi risali) e dopo un decrease_key (abbassi la priorità, potrebbe dover risalire)
void heap_sift_up (size_t i) {
    while(i > 0) {
        size_t parent = (i - 1) / 2;                            // calcola il padre dell'elemento in posizione i
        if(heap.arr[parent].dist <= heap.arr[i].dist) break;    // se il padre ha dist <= del figlio -> stop, invariante soddisfatta
        heap_swap(i, parent);                                   // altrimenti scambia figlio-padre
        i = parent;                                             // e continua a risalire con i = parent
    }
}

// HEAP_SIFT_DOWN: climbs down the heap
// Si usa dopo un pop (metti l'ultimo elemento in radice e poi riscendi confrontando root coi figli)
void heap_sift_down (size_t i) {
    while(1) {                                  // scende l'elemento in posizione i finché è minore dei suoi figli
        size_t left = 2*i + 1;                  // calcola indice figlio sinistro
        size_t right = 2*i + 2;                 // calcola indice figlio destro
        size_t smallest = i;


        // trova il più piccolo tra i, left e right
        if(left < heap.size && heap.arr[left].dist < heap.arr[smallest].dist) {
            smallest = left;
        }
        if(right < heap.size && heap.arr[right].dist < heap.arr[smallest].dist) {
            smallest = right;
        }

        if(smallest == i) break;                // se i è già il più piccolo -> stop
        heap_swap(i, smallest);                 // altrimenti scambia con il figlio più piccolo tra left e right
        i = smallest;                           // e continua a scnedere con i = smallest
    }
}

// HEAP_PUSH: adds a new node in the heap
void heap_push (int32_t node, int32_t dist) {
    size_t i = heap.size++;                     // inserisce in coda (i = size) ed incrementa size
    heap.arr[i].node = node;                    // Scrive la coppia (node, dist)
    heap.arr[i].dist = dist;
    heap.pos[node] = (int32_t)i;                // e aggiorna pos[node]= i
    heap_sift_up(i);                            // ripristina l'invariante salendo
}

// HEAP_POP: pops the node with the minimum distance
Heap_node heap_pop (void) {
    Heap_node root = heap.arr[0];               // salva la radice (nodo con dist minimo)
    heap.pos[root.node] = -1;                   // segnala che root non è piu nella coda

    heap.size--;                                // rimuove "logicamente" l'ultimo elemento
    if(heap.size > 0) {                         // se dopo la rimozione l'heap non è vuoto, dobbiamo riempire il buco in cima e ristabilire l'ordine
        heap.arr[0] = heap.arr[heap.size];      // sposto l'ultimo elemento alla radice per non lasciare buchi
        heap.pos[heap.arr[0].node] = 0;         // aggiorno la posizione del nodo che abbiamo appena portato in cima
        heap_sift_down(0);                      // chiamiamo sift_down per sistemare il min heap
    }
    return root;    // ritorna il minimo estratto (nodo + distanza)
}

// DECREASE_KEY: updates the distance of a node which is already inside the heap
void heap_decrease_key (int32_t node, int32_t new_dist) {
    int32_t i = heap.pos[node];                 // trovo l'indice da aggiornare in O(1) grazie alla tabella pos[]
    if (i == -1) return;                        // node outside the heap (non c'è nulla da fare)

    heap.arr[i].dist = new_dist;                // aggiorna la chiave (dist), ovvero la priorità del nodo
    heap_sift_up((size_t)i);                    // risaliamo finche non è rispettata la proprietà del min heap
}

// HEAP_EMPTY: checks if the heap is empty
int heap_empty (void) {
    return heap.size == 0;      // se ritorna 1, vuol dire che heap.size = 0, l'heap è vuoto e quindi non c'è nulla da estrarre
}


void travel_cost (int sx_in, int sy_in, int dx_in, int dy_in) {
    // validazioni base
    if(!map || !air_head || !heap.arr || !sp_dist || !sp_done) { puts("-1"); return; }        // controllo che tutto lo stato necessario esista
    if(sx_in < 0 || sy_in < 0 || dx_in < 0 || dy_in < 0) { puts("-1"); return; }             // controllo che i parametri di ingresso siano maggiori di zero

    // cast a size_t di tutti i parametri in input
    size_t sx = (size_t)sx_in;
    size_t sy = (size_t)sy_in;
    size_t dx = (size_t)dx_in;
    size_t dy = (size_t)dy_in;

    if(!in_bounds(sx, sy) || !in_bounds(dx, dy)) { puts("-1"); return; }  // controllo dei bounds dei due esagoni

    if(sx == dx && sy == dy) { puts("0"); return; }     // caso banale: sorgente = destinazione

    // reset heap e buffer di lavoro
    heap_reset();                           // azzero heap riusabile

    size_t n_cells = n_cols * n_rows;
    for(size_t i = 0; i < n_cells; i++) {   // inizializzazione vettori di lavoro
        sp_dist[i] = INF32;                 // set up distanze a infinito
        sp_done[i] = 0;                     // nessun nodo finalizzato
    }

    // calcolo degli indici 1d di sorgente (s) e destinazione (d)
    int32_t s = (int32_t)idx(sx, sy);
    int32_t d = (int32_t)idx(dx, dy);


    /* --- DIJKSTRA --- */
    // set up Dijkstra: la sorgente parte a 0 ed entra nella coda
    sp_dist[s] = 0;
    heap_push(s, 0);

    // loop principale Dijkstra
    while(!heap_empty()) {
        Heap_node hn = heap_pop();  // estrai il nodo con distanza minima stimata
        int32_t u = hn.node;

        if(sp_done[u]) continue;    // se già finalizzato, salta (key update concorrenti)

        sp_done[u] = 1;             // ora u ha distanza definitiva
        if(u == d) break;           // trovato il costo ottimo per la destinazione

        // ricavo coordinate 2d (col, row) per enumerare i vicini via terra
        size_t uc = (size_t)(u % (int32_t)n_cols);
        size_t ur = (size_t)(u / (int32_t)n_cols);

        /* --- Espansione via Terra --- */
        cost_t ground = map[u].cost;                            // costo di uscita via terra da u
        if(ground > 0) {                                        // se il costo di u = 0, non puoi uscire via terra
            size_t ncols[6], nrows[6];
            int k = neighbors_odd_r(uc, ur, ncols, nrows);      // fino a 0..6 vicini validi
            for (int i = 0; i < k; ++i) {
                int32_t v = (int32_t)idx(ncols[i], nrows[i]);
                int32_t alt = hn.dist + (int32_t)ground;        // paga l'uscita da u
                if(alt < sp_dist[v]) {                          // rilassamento standard
                    sp_dist[v] = alt;
                    if(heap.pos[v] == -1) {
                        heap_push(v, alt);                      // nuovo in priority queue
                    } else {
                        heap_decrease_key(v, alt);              // migliora priorità
                    }

                }

            }

        /* --- Espansione via Aria --- */

            for(int32_t e = air_head[u]; e != NO_EDGE; e = air_edges[e].next) {
                int32_t v = air_edges[e].to;                        // destinazione dell'arco aereo
                int32_t w = (int32_t)air_edges[e].cost;             // costo della rotta aerea u->v
                int32_t alt = hn.dist + w;                          // paghi la rotta e non il terreno di u
                if(alt < sp_dist[v]) {                              // rilassamento
                    sp_dist[v] = alt;
                    if(heap.pos[v] == -1) {
                            heap_push(v, alt);
                        } else {
                            heap_decrease_key(v, alt);
                        }
                }
            }
        }
    }
    /* --- Risulato --- */
        // Il costo di uscita della destinazione è già automaticamente escluso, perché i pesi si pagano quando si ESCE da un nodo, non quando si entra
        if(sp_dist[d] == INF32) {
            puts("-1");                     // irraggiungibile
        } else {
            printf("%d\n", sp_dist[d]);     // stampa del costo minimo
        }

}
