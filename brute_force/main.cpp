#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <stack>
#define INF 9999
using namespace std;
int** edge_list;
int number_of_cities;

void karp() {
    // initialize the DP table
// best[visited][last] = the cheapest way to reach the state where:
// - "visited" encodes the set of visited vertices other than N-1
// - "last" is the number of the vertex visited last

    vector< vector<int> > best( 1<<(number_of_cities-1), vector<int>( number_of_cities, INT_MAX ) );

// use DP to solve all states
// note that whenever we solve a particular state,
// the smaller states we need have already been solved

    for (int visited = 1; visited < (1<<(number_of_cities-1)); ++visited) {
        for (int last = 0; last < (number_of_cities-1); ++last) {

            // last visited vertex must be one of the visited vertices
            if (!(visited & 1<<last)) continue;

            // try all possibilities for the previous vertex,
            // pick the best among them
            if (visited == 1 << last) {
                // previous vertex must have been N-1
                best[visited][last] = edge_list[number_of_cities-1][last];
            } else {
                // previous vertex was one of the other ones in "visited"
                int prev_visited = visited ^ 1<<last;
                for (int prev = 0; prev < number_of_cities-1; ++prev) {
                    if (!(prev_visited & 1<<prev)) continue;
                    best[visited][last] = min(
                            best[visited][last],
                            edge_list[last][prev] + best[prev_visited][prev]
                    );
                }
            }
        }
    }

// use the precomputed path lengths to choose the cheapest cycle
    int answer = INT_MAX;
    for (int last=0; last<number_of_cities-1; ++last) {
        answer = min(
                answer,
                edge_list[last][number_of_cities-1] + best[ (1<<(number_of_cities-1))-1 ][last]
        );
    }
}
void readFile(string filename) {
    int representation_type;
    int edge_value;
    fstream file_stream;
    file_stream.open(filename);
    if (file_stream.is_open()) {
        file_stream >> representation_type;
        if(representation_type == 1) {
            file_stream >> number_of_cities;
            edge_list = new int*[number_of_cities];
            for(int i = 0; i < number_of_cities; i++) {
                edge_list[i] = new int[number_of_cities];
            }
            for(int i = 0; i < number_of_cities; i++) {
                for(int j = 0; j < number_of_cities; j++) {
                   file_stream >> edge_list[i][j];
                }
            }
            for(int i = 0; i < number_of_cities; i++) {
                for(int j = 0; j < number_of_cities; j++) {
                    cout << edge_list[i][j] << " ";
                }
                cout << endl;
            }
        }

    }
    file_stream.close();
}

int brute_force_tsp(string file_name, int iteration_amount, int optimal_solution) {
    vector<int> path_to_be_returned;
    int value_of_path_to_be_returned;
    fstream file_stream_results;
    file_stream_results.open("wyniki.csv", ios_base::app);
    file_stream_results << file_name << " " << iteration_amount << " " << optimal_solution << endl;
    for(int i = 0; i < iteration_amount; i++) {
        auto begin = std::chrono::system_clock::now();
        vector<int> path;
        int result;
        int start = 0;
        vector<int> nodes;
        for (int i = 0; i < number_of_cities; i++) {
            if (i != start) nodes.push_back(i);
        }
        int n = nodes.size();
        int shortest_path = INT_MAX;
        while (1) {
            int path_sum = 0;
            int j = start;
            path.clear();
            path.push_back(start);
            for (int i = 0; i < n; i++) {
                path_sum += edge_list[j][nodes[i]];
                path.push_back(nodes[i]);
                j = nodes[i];
            }
            path_sum += edge_list[j][start];
            if(path_sum < shortest_path) path_to_be_returned = path;
            shortest_path = min(shortest_path, path_sum);
            value_of_path_to_be_returned = shortest_path;
            if (!(next_permutation(nodes.begin(), nodes.end()))) break;
        }
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        file_stream_results << elapsed << endl;
    }
    for (int i = 0; i < number_of_cities; ++i) {
        delete[] edge_list[i];
    }
    delete[] edge_list;
    for(int i = 0; i < path_to_be_returned.size(); i++) {
        file_stream_results << path_to_be_returned[i] << " ";
    }
    file_stream_results << "0" << " ";
    file_stream_results << endl;
    file_stream_results.close();
    return value_of_path_to_be_returned;
 }




int main() {
    int znak;
    int number_of_instances = 1;
    string* instance_names = new string[number_of_instances];
    int* iteration_amounts = new int[number_of_instances];
    int* optimal_solutions = new int[number_of_instances];
    fstream initstream;
    initstream.open("init.txt");
    int i = 0;
    while(!initstream.eof()) {
        initstream >> instance_names[i];
        initstream >> iteration_amounts[i];
        initstream >> optimal_solutions[i];
        i++;
    }
    cout << "------------------------------------------------------" << std::endl;
    cout << "TSP PROBLEM SOLVER 1.0 [BRUTE FORCE]\n";
    cout << "------------------------------------------------------\n";
    cout << "Loading init" << endl;

     for(int i = 0; i < number_of_instances; i++) {
         readFile(instance_names[i]);
          brute_force_tsp(instance_names[i], iteration_amounts[i], optimal_solutions[i]);
         //karp();
      }
    return 0;

}
