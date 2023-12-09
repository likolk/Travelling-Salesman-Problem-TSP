#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <random>
#include <chrono>
#include <array>
#include <cstdlib> 
#include <ctime>


// the 10 benchmarks
const char *files[10] = {
    "../cities/ch130.tsp",
    "../cities/d198.tsp",
    "../cities/eil76.tsp",
    "../cities/fl1577.tsp",
    "../cities/kroA100.tsp",
    "../cities/lin318.tsp",
    "../cities/pcb442.tsp",
    "../cities/pr439.tsp",
    "../cities/rat783.tsp",
    "../cities/u1060.tsp"
};


void parse_file(const char* file_name, std::vector<std::pair<double, double> >& coords, int& best_known_route){
    std::ifstream file(file_name, std::ifstream::in);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    for(int i = 0; i < 5; ++i) {
        std::string line;
        getline(file, line);
    }

    std::string line;
    getline(file, line);
    std::stringstream best_known_ss(line);
    int best_known_dist;
    for(int i = 0; i < 3; ++i) {
        getline(best_known_ss, line, ' ');
    }

    try {
        best_known_dist = std::stoi(line);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
        std::cerr << "Line causing the issue: " << line << std::endl;
        file.close();
        return;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << std::endl;
        std::cerr << "Line causing the issue: " << line << std::endl;
        file.close();
        return;
    }

    best_known_route = best_known_dist;

    getline(file, line);
    while(getline(file, line)) {
        std::stringstream line_split(line);
        std::string number;
        getline(line_split, number, ' ');
        if(number != "EOF") {
            std::pair<double, double> coord;
            double values[2];
            for(int i = 0; i < 2; ++i) {
                getline(line_split, number, ' ');
                values[i] = std::stof(number, nullptr);
            }
            coord = std::make_pair(values[0], values[1]);
            coords.push_back(coord);
        }
    }
    file.close();
}


// Euclidean distance formula as given in the pdf.
int distance(std::pair<double, double> p1, std::pair<double, double> p2){
    double x, y;
    x = p1.first - p2.first;
    y = p1.second - p2.second;
    return static_cast<int>(sqrt(x*x + y*y)); // converts the expression to int
}

int calculate_route(const std::vector<std::pair<double, double> >& coords, const std::vector<int>& order){
    int route = 0;
    for (size_t i = 0; i < coords.size() - 1; ++i) {
        route += distance(coords[order[i]], coords[order[i + 1]]);
    }
    route += distance(coords[order[coords.size() - 1]], coords[order[0]]);
    return route;
}
std::vector<int> nearest_neighbor(const std::vector<std::pair<double, double> >& coords) {
    std::vector<int> order;
    order.reserve(coords.size());

    std::vector<bool> visited(coords.size(), false);
    order.push_back(0); // Start from the first city
    visited[0] = true;

    for (size_t i = 1; i < coords.size(); ++i) {
        int last = order.back();
        int nearest = -1;
        double min_dist = std::numeric_limits<double>::max();

        for (size_t j = 0; j < coords.size(); ++j) {
            if (!visited[j]) {
                double dist = distance(coords[last], coords[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    nearest = j;
                }
            }
        }

        order.push_back(nearest);
        visited[nearest] = true;
    }

    return order;
}


std::vector<int> two_opt(const std::vector<std::pair<double, double> >& coords, const std::vector<int>& order) {
    std::vector<int> new_order = order;
    int size = static_cast<int>(new_order.size());

    bool improved = true;
    int current_dist = calculate_route(coords, new_order);

    while (improved) {
        improved = false;
        for (int i = 1; i < size - 1; ++i) {
            for (int j = i + 1; j < size; ++j) {
                if (j - i == 1) continue; // No need to reverse if adjacent cities

                int new_dist = current_dist -
                    distance(coords[new_order[i - 1]], coords[new_order[i]]) -
                    distance(coords[new_order[j]], coords[new_order[(j + 1) % size]]) +
                    distance(coords[new_order[i - 1]], coords[new_order[j]]) +
                    distance(coords[new_order[i]], coords[new_order[(j + 1) % size]]);

                if (new_dist < current_dist) {
                    std::reverse(new_order.begin() + i, new_order.begin() + j + 1);
                    current_dist = new_dist;
                    improved = true;
                }
            }
        }
    }
    return new_order;
}

std::vector<int> two_and_half_opt(const std::vector<std::pair<double, double> >& coords, const std::vector<int>& order) {
    std::vector<int> new_order = order;
    int size = static_cast<int>(new_order.size());

    bool improved = true;
    int current_dist = calculate_route(coords, new_order);

    while (improved) {
        improved = false;
        for (int i = 1; i < size - 2; ++i) {
            for (int j = i + 2; j < size; ++j) {
                if (j - i == 1) continue; // No need to reverse if adjacent cities

                int new_dist = current_dist -
                    distance(coords[new_order[i - 1]], coords[new_order[i]]) -
                    distance(coords[new_order[j]], coords[new_order[(j + 1) % size]]) +
                    distance(coords[new_order[i - 1]], coords[new_order[j]]) +
                    distance(coords[new_order[i]], coords[new_order[(j + 1) % size]]);

                if (new_dist < current_dist) {
                    std::reverse(new_order.begin() + i, new_order.begin() + j + 1);
                    current_dist = new_dist;
                    improved = true;
                    break; // Exit inner loop if improvement is made
                }
            }
            if (improved) break; // Exit outer loop if improvement is made
        }
    }
    return new_order;
}




std::vector<int> three_opt(const std::vector<std::pair<double, double> >& coords, const std::vector<int>& order, double max_execution_time) {
    std::vector<int> new_order = order;
    int size = static_cast<int>(new_order.size());

    bool improved = true;
    int iteration = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    while (improved) {
        improved = false;
        for (int i = 1; i < size - 3; ++i) {
            for (int j = i + 2; j < size - 1; ++j) {
                for (int k = j + 2; k < size; ++k) {
                    int d0 = distance(coords[new_order[i - 1]], coords[new_order[i]]);
                    int d1 = distance(coords[new_order[j - 1]], coords[new_order[j]]);
                    int d2 = distance(coords[new_order[k - 1]], coords[new_order[k]]);

                    int delta = 0;
                    delta += distance(coords[new_order[i - 1]], coords[new_order[j - 1]]) - d0;
                    delta += distance(coords[new_order[i]], coords[new_order[k - 1]]) - d1;
                    delta += distance(coords[new_order[j]], coords[new_order[k]]) - d2;

                    if (delta < 0) {
                        std::reverse(new_order.begin() + i, new_order.begin() + j);
                        std::reverse(new_order.begin() + j, new_order.begin() + k);
                        improved = true;
                    }
                }
            }
        }

        ++iteration;

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        if (elapsed_time.count() > max_execution_time) {
            std::cout << "Maximum execution time reached. Stopping the algorithm." << std::endl;
            break;
        }
    }

    return new_order;
}

std::vector<int> perturb_solution(const std::vector<int>& order) {
    std::vector<int> perturbed_order = order;
    // Implement a perturbation mechanism (perturbing subsets or fragments)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> index_dist(0, perturbed_order.size() - 1);
    int start_index = index_dist(gen);
    int end_index = index_dist(gen);

    if (start_index > end_index) {
        std::swap(start_index, end_index);
    }

    // Perturb the subset
    std::reverse(perturbed_order.begin() + start_index, perturbed_order.begin() + end_index + 1);

    return perturbed_order;
}



// Basic Iterated Local Search (ILS)
std::vector<int> iterated_local_search(const std::vector<std::pair<double, double> >& coords, const std::vector<int>& order) {
    std::vector<int> current_order = order;
    int current_distance = calculate_route(coords, current_order);
    const int MAX_ITERATIONS = 30;

    for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
        std::vector<int> perturbed_order = perturb_solution(current_order);
        perturbed_order = two_and_half_opt(coords, perturbed_order);
        int perturbed_distance = calculate_route(coords, perturbed_order);

        // Accept the perturbed solution if it's better than the current one
        if (perturbed_distance < current_distance) {
            current_order = perturbed_order;
            current_distance = perturbed_distance;
        }
    }
    return current_order;
}

int main() {
    std::srand(100);
    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < 10; ++i) {
        int best_known_route;
        std::vector<std::pair<double, double> > coords;
        parse_file(files[i], coords, best_known_route);

        std::vector<int> best_order = nearest_neighbor(coords);

        // Perform Iterated Local Search (2.5-opt)
        best_order = iterated_local_search(coords, best_order);
        int min_route = calculate_route(coords, best_order);

        std::cout << files[i] << " 2.5-opt: " << min_route << " Best known: " << best_known_route << std::endl;
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}
