#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <stack>
#include <time.h>
#include <ctime>
#include <random>
#include<utility>

#define INF 9999
using namespace std;
int** edge_list;
int** stateArray;
int** getNodeArray;
int number_of_cities;


int calcPath(vector<int> path) {
    int cost = 0;

    for (int i = 0; i < path.size() - 1; ++i)
    {
        cost += edge_list[path[i]][path[i + 1]];
    }
    cost += edge_list[path[number_of_cities - 1]][path[0]];

    return cost;
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

int calculateCNN(){

    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_int_distribution<int> d(1, number_of_cities-1);

    vector<int> path;
    for( int k=0; k < number_of_cities; k++)
        path.push_back(k);
    path.push_back(0);

    int sum = 0;
    int N = 10000;
    for(int i = 0; i < N; i++){
        std::shuffle(path.begin()+1 ,path.end()-1, g);
        sum += calcPath(path);
    }
    sum = (double) sum / (double) N;
    return sum;
}

void ACO(string filename, int iteration_amount) {
    for(int z = 0; z < iteration_amount; z++) {
        auto begin = std::chrono::system_clock::now();
        fstream file_stream;
        file_stream.open("wyniki.csv", ios_base::app);
        file_stream << filename << endl;
        std::random_device rd;
        std::mt19937 g(rd());
        std::uniform_real_distribution<double> num(0.0, 1.0);
        std::uniform_int_distribution<int> node(0, number_of_cities - 1);
        vector<int> best_path;
        int best_cost = -1;

        double alpha = 1.0;
        double beta = 2.5;
        double evaporation = 0.5;

        int Cnn = calculateCNN();
        double init_pher = (double) number_of_cities / (double) Cnn;
        int iterations = 100;
        double pher = 100.0;

        vector<vector<double>> pheromones;
        for (int i = 0; i < number_of_cities; i++) {
            vector<double> ph;
            for (int k = 0; k < number_of_cities; k++)
                ph.push_back(init_pher);
            pheromones.push_back(ph);
        }

        vector<vector<double>> pheromones_memory;
        for (int i = 0; i < number_of_cities; i++) {
            vector<double> ph;
            for (int k = 0; k < number_of_cities; k++)
                ph.push_back(0.0);
            pheromones_memory.push_back(ph);
        }
        int i = 0;
        while (i++ < iterations) {
            std::cout << "ITERACJA: " << i << std::endl;

            for (int m = 0; m < number_of_cities; m++) {

                int start_node = m;
                vector<int> path;
                path.push_back(start_node);
                vector<int> left_vertices;

                for (int y = 0; y < number_of_cities; y++) {
                    if (y != start_node)
                        left_vertices.push_back(y);
                }

                for (int x = 0; x < number_of_cities - 1; x++) {
                    int cn = path.back();

                    double sum_factor = 0.0;
                    for (vector<int>::iterator it = left_vertices.begin(); it != left_vertices.end(); it++)
                        sum_factor += (pow(pheromones[cn][*it], alpha)) *
                                      (pow((double) 1.0 / (double) edge_list[cn][*it], beta));

                    vector<double> probabil;

                    for (vector<int>::iterator ir = left_vertices.begin(); ir != left_vertices.end(); ir++) {
                        double p = ((pow(pheromones[cn][*ir], alpha)) *
                                    (pow((double) 1.0 / (double) edge_list[cn][*ir], beta))) / sum_factor;
                        probabil.push_back(p);
                    }

                    partial_sum(probabil.begin(), probabil.end(), probabil.begin());

                    double p = num(g);

                    int k = 0;

                    for (vector<double>::iterator ix = probabil.begin(); ix != probabil.end(); ix++) {

                        if (p < *ix) {

                            int wybrany = left_vertices[k];
                            path.push_back(wybrany);
                            left_vertices.erase(left_vertices.begin() + k);

                            break;
                        }
                        k++;
                    }

                }


                int L = calcPath(path);
                if (L < best_cost || best_cost == -1) {
                    best_cost = L;
                    best_path = path;
                    i = 0;
                }

                double ph = pher / (double) L;

                for (int r = 0; r < path.size() - 1; r++)
                    pheromones_memory[path[r]][path[r + 1]] += ph;

            }

            for (int i = 0; i < number_of_cities; i++) {
                for (int k = 0; k < number_of_cities; k++) {
                    pheromones[i][k] = (pheromones[i][k] * evaporation) + pheromones_memory[i][k];
                    pheromones_memory[i][k] = 0.0;
                }

            }

        }
        for (int z = 0; z < best_path.size(); z++) {
            std::cout << best_path[z] << " ";
        }
        cout << best_cost;
        std::cout << std::endl;
        best_path.pop_back();
        vector<int>::iterator itt = find(best_path.begin(), best_path.end(), 0);
        std::rotate(best_path.begin(), itt, best_path.end());
        best_path.push_back(0);
        vector<int> solution_and_cost(best_path);
        solution_and_cost.push_back(best_cost);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        file_stream << elapsed << "<- czas" << " " << best_cost << "<- koszt" << " CAS <-- typ" << endl;
        file_stream.close();
    }
}

void ACO2(string filename, int iteration_amount) {
    for(int z = 0; z < iteration_amount; z++) {
        auto begin = std::chrono::system_clock::now();


        fstream file_stream;
        file_stream.open("wyniki.csv", ios_base::app);
        file_stream << filename << endl;
        std::random_device rd;
        std::mt19937 g(rd());
        std::uniform_real_distribution<double> num(0.0, 1.0);
        std::uniform_int_distribution<int> node(0, number_of_cities - 1);
        vector<int> best_path;
        int best_cost = -1;

        double alpha = 0.0;
        double beta = 2.5;
        double evaporation = 0.5; // ro

        int Cnn = calculateCNN();
        double init_pher = (double) number_of_cities / (double) Cnn;
        int iterations = 100;
        double pher = 100.0;

        vector<vector<double>> pheromones;
        for (int i = 0; i < number_of_cities; i++) {
            vector<double> ph;
            for (int k = 0; k < number_of_cities; k++)
                ph.push_back(init_pher);
            pheromones.push_back(ph);
        }

        vector<vector<double>> pheromones_memory;
        for (int i = 0; i < number_of_cities; i++) {
            vector<double> ph;
            for (int k = 0; k < number_of_cities; k++)
                ph.push_back(0.0);
            pheromones_memory.push_back(ph);
        }
        int i = 0;
        while (i++ < iterations) {
            std::cout << "ITERACJA: " << i << std::endl;

            for (int m = 0; m < number_of_cities; m++) {

                int start_node = m;
                vector<int> path;
                path.push_back(start_node);
                vector<int> left_vertices;

                for (int y = 0; y < number_of_cities; y++) {
                    if (y != start_node)
                        left_vertices.push_back(y);
                }
                // szukanie sciezki
                for (int x = 0; x < number_of_cities - 1; x++) {
                    int cn = path.back();
                    double sum_factor = 0.0;
                    for (vector<int>::iterator it = left_vertices.begin(); it != left_vertices.end(); it++)
                        sum_factor += (pow(pheromones[cn][*it], alpha)) *
                                      (pow((double) 1.0 / (double) edge_list[cn][*it], beta));

                    vector<double> probabil;

                    for (vector<int>::iterator ir = left_vertices.begin(); ir != left_vertices.end(); ir++) {
                        double p = ((pow(pheromones[cn][*ir], alpha)) *
                                    (pow((double) 1.0 / (double) edge_list[cn][*ir], beta))) / sum_factor;
                        probabil.push_back(p);
                    }

                    partial_sum(probabil.begin(), probabil.end(), probabil.begin());

                    double p = num(g);

                    int k = 0;

                    for (vector<double>::iterator ix = probabil.begin(); ix != probabil.end(); ix++) {

                        if (p < *ix) {

                            int wybrany = left_vertices[k];
                            path.push_back(wybrany);
                            left_vertices.erase(left_vertices.begin() + k);

                            break;
                        }
                        k++;
                    }

                }


                int L = calcPath(path);
                if (L < best_cost || best_cost == -1) {
                    best_cost = L;
                    best_path = path;
                    i = 0;
                }

                double ph = pher;

                for (int r = 0; r < path.size() - 1; r++)
                    pheromones_memory[path[r]][path[r + 1]] += ph;

            }

            for (int i = 0; i < number_of_cities; i++) {
                for (int k = 0; k < number_of_cities; k++) {
                    pheromones[i][k] = (pheromones[i][k] * evaporation) + (pheromones_memory[i][k]/edge_list[i][k]);
                    pheromones_memory[i][k] = 0.0;
                }

            }

        }
        for (int z = 0; z < best_path.size(); z++) {
            std::cout << best_path[z] << " ";
        }
        cout << best_cost;
        std::cout << std::endl;
        best_path.pop_back();
        vector<int>::iterator itt = find(best_path.begin(), best_path.end(), 0);
        std::rotate(best_path.begin(), itt, best_path.end());
        best_path.push_back(0);
        vector<int> solution_and_cost(best_path);
        solution_and_cost.push_back(best_cost);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        file_stream << elapsed << "<- czas" << " " << best_cost << "<- koszt" << " QAS <-- typ" << endl;
        file_stream.close();
    }
}

void ACO3(string filename, int iteration_amount) {
    for(int z = 0; z < iteration_amount; z++) {
        auto begin = std::chrono::system_clock::now();


        fstream file_stream;
        file_stream.open("wyniki.csv", ios_base::app);
        file_stream << filename << endl;
        std::random_device rd;
        std::mt19937 g(rd());
        std::uniform_real_distribution<double> num(0.0, 1.0);
        std::uniform_int_distribution<int> node(0, number_of_cities - 1);
        vector<int> best_path;
        int best_cost = -1;

        double alpha = 1.0;
        double beta = 2.5;
        double evaporation = 0.5; // ro

        int Cnn = calculateCNN();
        double init_pher = (double) number_of_cities / (double) Cnn;
        int iterations = 100;
        double pher = 100.0;

        vector<vector<double>> pheromones;
        for (int i = 0; i < number_of_cities; i++) {
            vector<double> ph;
            for (int k = 0; k < number_of_cities; k++)
                ph.push_back(init_pher);
            pheromones.push_back(ph);
        }

        vector<vector<double>> pheromones_memory;
        for (int i = 0; i < number_of_cities; i++) {
            vector<double> ph;
            for (int k = 0; k < number_of_cities; k++)
                ph.push_back(0.0);
            pheromones_memory.push_back(ph);
        }
        int i = 0;
        while (i++ < iterations) {
            std::cout << "ITERACJA: " << i << std::endl;

            for (int m = 0; m < number_of_cities; m++) {

                int start_node = m;
                vector<int> path;
                path.push_back(start_node);
                vector<int> left_vertices;

                for (int y = 0; y < number_of_cities; y++) { //
                    if (y != start_node)
                        left_vertices.push_back(y);
                }
                for (int x = 0; x < number_of_cities - 1; x++) {
                    int cn = path.back();

                    double sum_factor = 0.0;
                    for (vector<int>::iterator it = left_vertices.begin(); it != left_vertices.end(); it++)
                        sum_factor += (pow(pheromones[cn][*it], alpha)) *
                                      (pow((double) 1.0 / (double) edge_list[cn][*it], beta));

                    vector<double> probabil;

                    for (vector<int>::iterator ir = left_vertices.begin(); ir != left_vertices.end(); ir++) {
                        double p = ((pow(pheromones[cn][*ir], alpha)) *
                                    (pow((double) 1.0 / (double) edge_list[cn][*ir], beta))) / sum_factor;
                        probabil.push_back(p);
                    }

                    partial_sum(probabil.begin(), probabil.end(), probabil.begin());

                    double p = num(g);

                    int k = 0;

                    for (vector<double>::iterator ix = probabil.begin(); ix != probabil.end(); ix++) {

                        if (p < *ix) {

                            int wybrany = left_vertices[k];
                            path.push_back(wybrany);
                            left_vertices.erase(left_vertices.begin() + k);

                            break;
                        }
                        k++;
                    }

                }

                int L = calcPath(path);
                if (L < best_cost || best_cost == -1) {
                    best_cost = L;
                    best_path = path;
                    i = 0;
                }

                double ph = pher;

                for (int r = 0; r < path.size() - 1; r++)
                    pheromones_memory[path[r]][path[r + 1]] += ph;

            }

            for (int i = 0; i < number_of_cities; i++) {
                for (int k = 0; k < number_of_cities; k++) {
                    pheromones[i][k] = (pheromones[i][k] * evaporation) + (pheromones_memory[i][k]);
                    pheromones_memory[i][k] = 0.0;
                }

            }

        }
        for (int z = 0; z < best_path.size(); z++) {
            std::cout << best_path[z] << " ";
        }
        cout << best_cost;
        std::cout << std::endl;
        best_path.pop_back();
        vector<int>::iterator itt = find(best_path.begin(), best_path.end(), 0);
        std::rotate(best_path.begin(), itt, best_path.end());
        best_path.push_back(0);
        vector<int> solution_and_cost(best_path);
        solution_and_cost.push_back(best_cost);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        file_stream << elapsed << "<- czas" << " " << best_cost << "<- koszt" << " DAS <-- typ" << endl;
        file_stream.close();
    }
}

int main() {
    srand(time(NULL));

    int znak;
    int number_of_instances = 14;
    string* instance_names = new string[number_of_instances];
    int* iteration_amounts = new int[number_of_instances];
    fstream initstream;
    initstream.open("init.txt");
    int i = 0;
    while(!initstream.eof()) {
        initstream >> instance_names[i];
        initstream >> iteration_amounts[i];
        i++;
    }
    cout << "------------------------------------------------------" << std::endl;
    cout << "TSP PROBLEM SOLVER 3.0 [ACO]\n";
    cout << "------------------------------------------------------\n";
    cout << "Loading init" << endl;
    for(int i = 0; i < number_of_instances; i++) {
        readFile(instance_names[i]);
        ACO3(instance_names[i], iteration_amounts[i]);
    }
    return 0;

}

