#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>

// the 10 benchmarks
const char *files[10] = {
    "./data/ch130.tsp",
    "./data/d198.tsp",
    "./data/eil76.tsp",
    "./data/fl1577.tsp",
    "./data/kroA100.tsp",
    "./data/lin318.tsp",
    "./data/pcb442.tsp",
    "./data/pr439.tsp",
    "./data/rat783.tsp",
    "./data/u1060.tsp"
};


//get coordinates (coords) from file.   // a vector of pairs, taking two "double" type values
void parse_file(const char* file_name, std::vector<std::pair<double, double> > &coords, int &best_known_route){
    std::ifstream file(file_name, std::ifstream::in); // ifstream which takes a char pointer to a file_name buffer. std::ifstream::in is default from the documentation
    // get best_known_route
    for(int i=0;i<5;i++){
        std::string line;
        getline(file, line);
    }
    std::string line;
    getline(file, line);  // reads characters from file and places them into the "line" string.
    std::stringstream best_known_ss(line); // operates on strings.
    int best_known_dist;
    //parse the "best known line". Split line by spaces (' ')
    for(int i=0;i<3;i++){
        getline(best_known_ss, line, ' ');
    }
    best_known_dist = std::stoi(line); // parse line (basically converts "line" to an integer)
    best_known_route = best_known_dist;
    
    getline(file, line); // reads characters from file and places them into the "line" string.
    //get coords;
    while(getline(file, line)){ // while the file has not NULL number of characters.
        //parse coord line
        std::stringstream line_split(line); // operates on strings.
        std::string number;
        getline(line_split, number, ' ');
        if(number!="EOF"){//check if end of file
            //std::cout<<number<<std::endl;
            // use the pair class template to store the two (heterogeneous) objects, in this case, coordinates,
            // in a single unit, named coord.
            std::pair<double, double> coord;
            double values[2];
            for(int i=0;i<2;i++){
                getline(line_split, number, ' ');
                values[i] = std::stof(number, nullptr); // use stringtofloat stof to parse the string
                // "number" and as a float type
            } // take coord variable and assign to it the 2 values of the array, converted to a 
            // pair of x and y coordinates. Using make_pair function.
            coord = std::make_pair(values[0], values[1]); 
            coords.push_back(coord); // add it at the end of the vector
        }
    }
    // Here we use the close() function, to close the opened file named (apparently) file.
    file.close();
}
// Euclidean distance formula as given in the pdf.
int distance(std::pair<double, double> p1, std::pair<double, double> p2){
    double x, y;
    x = p1.first - p2.first;
    y = p1.second - p2.second;
    return static_cast<int>(sqrt(x*x + y*y)); // converts the expression to int
}


//idea of algorithm
//from starting point code goes to nearest point. Add distance to route. 
//then remove visited point from copy of coords vector. 
// return route             // a vector of pairs, taking two "double" type values
int min_route(std::vector<std::pair<double, double> > coords, int start_id){
    int route = 0, coords_size = coords.size();
    // a vector of pairs, taking two "double" type values
    std::pair<double, double> start_point = coords[start_id];
    // a vector of pairs, taking two "double" type values
    std::pair<double, double> curr_point = start_point;
    coords.erase(coords.begin()+start_id);
    for(int i=1; i<coords_size;i++){
        int min_dist = distance(curr_point, coords[0]); // calculate the euclidian distance between the current city we are right now and the coordinates at pos 0 and store it to the min_dist variable 
        int j = 0, id=0;
        for(auto c: coords){
            int dist = distance(c, curr_point);
            if(dist < min_dist){ // if the distance on the current city we are right now is greater than distance from position c to currect city
                min_dist = dist;
                id = j;
            }
            j += 1; 
        }
        route += min_dist;
        curr_point = coords[id];
        //remove visited courd
        coords.erase(coords.begin()+id);
    }
    route += distance(start_point, curr_point);
    return route;
}

int main(){
    auto start = std::chrono::steady_clock::now();
    // coords, that is a vector of pairs, taking two "double" type values
    std::vector<std::pair<double, double> > coords;
    // variation of nearest neighbour. Instead of checking nearest neighbour from only one starting point. Check every other possible starting point.
    // Code tweaks from 57.0 at the best case to 72 seconds every time it is run.
    for(int i=0;i<10;i++){
        int best_known_route;
        parse_file(files[i], coords, best_known_route); // parse file files[i]
        int min_r = min_route(coords, 0);
        for(int j=1;j<coords.size();j++){
            int r = min_route(coords, j);
            if(r < min_r){
                min_r = r;
            }
        }
        std::cout<<files[i]<<" route: "<<min_r<<" Best known: "<<best_known_route<<std::endl;
        coords.clear(); // destroys the coords vector content
    } // to return the total time needed to "iterate" between each vertex (in this case each city)
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}





/*
Reference Area:
https://www.baeldung.com/cs/tsp-dynamic-programming (for understanding mainly what TSP is and how it works, what its pseudocode looks like etc)
https://en.cppreference.com/w/cpp/chrono/steady_clock/now
*/

// K-Nearest Neighbour algorithm used to solve TSP

