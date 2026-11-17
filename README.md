# Travelling Salesman Problem

Heuristic solvers for the Travelling Salesman Problem in C++, written for the Artificial Intelligence course at Università della Svizzera italiana (USI).

## What is inside

| File | Approach |
| --- | --- |
| `TSP.cpp` | Nearest-neighbour construction |
| `solution.cpp` | Local search with 2-opt and 3-opt moves |
| `sixpercent.cpp` | Iterated local search |
| `Cities/` | TSPLIB instances (`eil76`, `kroA100`, `ch130`, `d198`, `lin318`, `pr439`, `pcb442`, `rat783`, `u1060`, `fl1577`) and two optimal tours |
| `Results.xls` | Tour lengths obtained on the instances |

## Build and run

```bash
g++ -O2 -std=c++17 solution.cpp -o tsp
./tsp
```

Each source file is a standalone program. Check the top of the file for how the instance is selected.
