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

vector<int> random_permutation(int _size)
{
    std::vector<int> temp;
    temp.reserve(_size);

    for (int i = 0; i < _size; i++)
    {
        temp.push_back(i);
    }

    random_shuffle(temp.begin(), temp.end());

    return temp;
}


int calculatePath(vector<int> path) {
    int cost = 0;

    for (int i = 0; i < path.size() - 1; ++i)
    {
        cost += edge_list[path[i]][path[i + 1]];
    }
    cost += edge_list[path[number_of_cities - 1]][path[0]];

    return cost;
}

double getProbability(int diff, double temperature) {
    return exp(diff / temperature);
}

double calculateTemperature() {

    vector<int> origin;

    int firstToSwap;
    int secondToSwap;
    int originCost = 0;
    int neighbourCost = 0;
    int delta = 0;
    int buffer = 0;

    for (int i = 0; i < 10000; i++)
    {
        do
        {
            firstToSwap = rand() % number_of_cities;
            secondToSwap = rand() % number_of_cities;
        } while (firstToSwap == secondToSwap);

        origin = random_permutation(number_of_cities);
        vector<int> neighbour(origin);

        std::swap(neighbour[firstToSwap], neighbour[secondToSwap]);

        delta = fabs(calculatePath(origin) - calculatePath(neighbour));
        buffer += delta;

    }

    buffer /= 10000;

    return (-1*buffer)/log(0.99);
}

void SA(string file_name, int iteration_amount, int optimal_solution, double timeF, double rate) {
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    fstream file_stream_results;
    file_stream_results.open("wyniki.csv", ios_base::app);
    vector<int> best;
    vector<int> permutation = random_permutation(number_of_cities);
    vector<int> next(permutation);
    std::clock_t start;
    double temperatureBuffer = calculateTemperature();
    vector<int> solution;
    int firstToSwap;
    int secondToSwap;
    double coolingRate = rate;
    double timeBound = timeF;
    double temperature =  temperatureBuffer;
   // double temperature =  1;
    cout << calculateTemperature() << endl;
    int result = 1 << 30;
    double time = 0;
    double foundTime = 0;
    start = std::clock();

    while(true)
    {
        while (temperature >= 0.1)
        {
            int steps =  3*number_of_cities;
            next = permutation;
            //cout << temperature << endl;
            int nextCost = calculatePath(next);
            for (int i = steps; i > 0; i--)
             {
                do
                {
                    firstToSwap = rand() % number_of_cities;
                    secondToSwap = rand() % number_of_cities;
                } while (firstToSwap == secondToSwap);

                std::swap(next[firstToSwap],next[secondToSwap]);

                nextCost = calculatePath(next);

                int difference = result - nextCost;
                 cout << getProbability(difference,temperature) << endl;
                if (difference > 0)
                {
                    result = nextCost;
                    best = next;
                    foundTime = (std::clock() - start) / (double)CLOCKS_PER_SEC;
                }


                if (difference > 0 || (difference < 0 && getProbability(difference,temperature) > dist(e2)))
                {
                    permutation = next;
                    break;
                }
                else
                {
                    std::swap(next[firstToSwap], next[secondToSwap]);
                 }

                time = (std::clock() - start) / (double)CLOCKS_PER_SEC;
                //cout << time << endl;
                if (time >= timeBound)
                {
                    file_stream_results << file_name << endl;
                    file_stream_results << timeF << endl;
                    cout << "Droga: ";
                    for (int d = 0; d < number_of_cities; d++)
                    {
                        file_stream_results << best[d] << " ";
                        cout << best[d] << " ";
                    }
                    file_stream_results << endl;
                    cout << "\nKoszt: " << result << endl;
                    file_stream_results << result << endl;
                    cout << "Znaleziono po: " << foundTime << " s " << endl;
                    file_stream_results << foundTime << endl;
                    cout << "Temperatura koncowa: "<< temperature << endl;
                    file_stream_results << temperature << endl;
                    cout << "Temperatura poczatkowa: "<< temperatureBuffer << endl;
                    cout << endl;
                    file_stream_results.close();
                    return;
                }
                }
            temperature *= coolingRate;
            }

        temperature = temperatureBuffer;
       // temperature = 1;
        permutation = random_permutation(number_of_cities);

    }




}

void karp(string file_name, int iteration_amount, int optimal_solution) {
    int valueToBeWritten;
    bool firstIteration = true;
    stack<int> path;
    fstream file_stream_results;
    file_stream_results.open("wyniki.csv", ios_base::app);
    file_stream_results << file_name << " " << iteration_amount << " " << optimal_solution << endl;
    for(int i = 0; i < iteration_amount; i++) {
        if(i == 1) firstIteration = false;
        auto begin = std::chrono::system_clock::now();
        int count = number_of_cities;
        stateArray = new int *[1 << count];
        getNodeArray = new int *[1 << count];
        int temp, result, lastVertex = 0;

        for (int i = 0; i < (1 << count); i++) {
            stateArray[i] = new int[count];
            getNodeArray[i] = new int[count];
        }
        for (int i = 0; i < (1 << count); i++) {
            for (int j = 0; j < count; j++) {
                stateArray[i][j] = INF;
                getNodeArray[i][j] = INF;
            }
        }

        for (int i = 0; i < count; i++) {
            stateArray[(1 << i) | 1][i] = edge_list[0][i];
            getNodeArray[(1 << i) | 1][i] = 0;
        }

        for (int bitMask = 0; bitMask < 1 << count; bitMask++) {
            for (int v = 0; v < count; v++) {
                if (!(bitMask & (1 << v)))
                    continue;

                for (int j = 0; j < count; j++) {
                    if (!(bitMask & (1 << j))) {
                        if (stateArray[bitMask][v] + edge_list[v][j] <= stateArray[bitMask | (1 << j)][j]) {
                            getNodeArray[bitMask | (1 << j)][j] = v;
                            stateArray[bitMask | (1 << j)][j] = stateArray[bitMask][v] + edge_list[v][j];
                        }
                    }
                }
            }
        }

        result = INF;
        temp = 0;

        for (int i = 0; i < count; i++) {
            temp = stateArray[(1 << count) - 1][i] + edge_list[i][0];
            if (temp <= result) {
                result = temp;
                lastVertex = i;
            }
        }

        valueToBeWritten = result;

        if(firstIteration) {
            int bitMask = (1 << count) - 1;
            int prev;

            path.push(lastVertex);
            for (int i = count - 2; i > 0; i--) {
                prev = getNodeArray[bitMask][lastVertex];
                path.push(prev);
                bitMask = bitMask ^ (1 << lastVertex);
                lastVertex = prev;
            }
        }

        cout << endl << result;
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        file_stream_results << elapsed << endl;
        for (int i = 0; i < (1 << number_of_cities); ++i) {
            delete[] stateArray[i];
        }
        delete[] stateArray;

        for (int i = 0; i < (1 << number_of_cities); ++i) {
            delete[] getNodeArray[i];
        }
        delete[] getNodeArray;
    }

    for (int i = 0; i < number_of_cities; ++i) {
        delete[] edge_list[i];
    }
    delete[] edge_list;
    int* path_holder = new int[number_of_cities];
    int i = 0;
    while(path.size()) {
        path_holder[i] = path.top();
        path.pop();
        i++;
    }
    file_stream_results << "0" << " ";
    path_holder[number_of_cities-1] = 0;
    for(int i = 0; i < number_of_cities; i++) {
        file_stream_results << path_holder[i] << " ";
    }

    file_stream_results << endl;
    file_stream_results << valueToBeWritten << endl;

    delete[] path_holder;

    file_stream_results << endl;
    file_stream_results.close();

}

std::vector<int> ShuffleVertices(std::vector<int> vector)
{
    std::vector<int> temp = vector;

    int x, y;

    do
    {
        x = rand() % number_of_cities;
        y = rand() % number_of_cities;
    }
    while (x == y);

    std::swap(temp.at(x), temp.at(y));

    return temp;
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

int main() {
    srand(time(NULL));

    int znak;
    int number_of_instances = 1;
    string* instance_names = new string[number_of_instances];
    int* iteration_amounts = new int[number_of_instances];
    int* optimal_solutions = new int[number_of_instances];
    int* times = new int[number_of_instances];
    fstream initstream;
    initstream.open("init.txt");
    int i = 0;
    while(!initstream.eof()) {
        initstream >> instance_names[i];
        initstream >> iteration_amounts[i];
        initstream >> optimal_solutions[i];
        initstream >> times[i];
        i++;
    }
    cout << "------------------------------------------------------" << std::endl;
    cout << "TSP PROBLEM SOLVER 4.0 [SA/KARP]\n";
    cout << "------------------------------------------------------\n";
    cout << "Loading init" << endl;

     /*for(int i = 0; i < number_of_instances; i++) {
         readFile(instance_names[i]);
         karp(instance_names[i], iteration_amounts[i], optimal_solutions[i]);
      }*/
    for(int i = 0; i < number_of_instances; i++) {
        readFile(instance_names[i]);
        SA(instance_names[i], iteration_amounts[i], optimal_solutions[i], times[i], 0.99);
    }
    return 0;

}
