# MPI-Smart-Arrays

Proof of concept thingy for writing MPI parallel programs without needing to worry about the MPI stuff. Probably buggy, just a proof of concept.

## Conways Game of Life:
```c++
#include "mpi-smart-arrays.hpp"

#include <cstdlib>

using cell = int;

struct Initializer {
    void operator()(cell* x, std::size_t global_idx, std::size_t global_size,
                    std::size_t elem_size) {
        for (size_t i = 0; i < elem_size; i++) {
            x[i] = rand() % 2;
        }
    }
};

struct Step {
    void operator()(cell* x, const cell* window, std::size_t window_size,
                    std::size_t global_idx, std::size_t global_size,
                    std::size_t elem_size) {
        const cell* l_cells = window - elem_size;
        const cell* r_cells = window + elem_size;
        const int ny = static_cast<int>(elem_size);
        for (int y = 0; y < ny; y++) {
            const int y_up = (y + 1) % ny;
            const int y_down = ((y - 1) + ny) % ny;
            const int nalive = x[y_up] + x[y_down] + l_cells[y_up] +
                               l_cells[y] + l_cells[y_down] + r_cells[y_up] +
                               r_cells[y] + r_cells[y_down];
            x[y] =
                (x[y] == 1) ? ((nalive >= 2) && (nalive <= 3)) : (nalive == 3);
        }
    }
};

int main() {
    MPI_Init(NULL, NULL);

    int seed = 3082002;
    srand(seed);

    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    std::size_t grid_x = 50;
    std::size_t grid_y = 50;

    smartarray::SmartArray<cell> arr(grid_x, grid_y, MPI_COMM_WORLD);
    arr.transform(Initializer());

    arr.print();

    for (int step = 0; step < 100; step++)
        arr.transform_window(Step(), 1);

    if (!comm_rank)
        std::cout << std::endl;
    arr.print();

    MPI_Finalize();
    return 0;
}
```