#include "Data.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <chrono> 

using namespace std;

struct solution{
    vector<int> sequence;
    double cost = 0.0;
};

struct subsequence{
    double T, C;
    int W, first, last; 
};

solution construcao(Data& data, double alpha){
    solution s;
    vector<int> CL;

    struct Candidato{
            int no;
            double distancia;
    };

    vector<Candidato>candidatos;
    vector<Candidato>RCL;

    s.sequence.push_back(1);

    for(int i = 2; i <= data.getDimension(); i++){
        CL.push_back(i);
    }

    while(!CL.empty()){
        candidatos.clear();
        RCL.clear();

        int r = s.sequence.back(); // no atual
        
        for(int k : CL){
            Candidato candidato;
            candidato.no = k;
            candidato.distancia = data.getDistance(r, k);
            candidatos.push_back(candidato);
        }

        sort(candidatos.begin(), candidatos.end(),
            [](const Candidato& a, const Candidato& b){
                return a.distancia < b.distancia;
            });

        int quantidadeRCL = (int)ceil(alpha * candidatos.size());
        
        for(int i = 0; i < quantidadeRCL; i++){
            RCL.push_back(candidatos[i]);
        }

        int indiceEscolhido = rand() % RCL.size();
        Candidato c = RCL[indiceEscolhido];

        s.sequence.push_back(c.no);

        CL.erase(find(CL.begin(), CL.end(), c.no));
    }

    s.sequence.push_back(1);
    return s;
}

int main(int argc, char** argv){
    auto start = chrono::high_resolution_clock::now();

    srand(time(NULL));

    auto data = Data(argc, argv[1]);
    data.read();
    size_t n = data.getDimension();

    cout << "Dimension: " << n << endl;

    double alpha = (double)rand() / RAND_MAX;

    solution s = construcao(data, alpha);

    cout << "\nExecutando Construcao para " << n << " cidades..." << endl;
    
    cout << "Melhor solucao encontrada: 0 -> ";
    for (int i = 0; i < (int)n; i++){
        cout << s.sequence[i] << " -> ";
    }
    cout << s.sequence[0] << " -> 0" << endl;

    cout << "Custo da melhor solucao: " << s.cost << endl;

    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration_sec = end - start;
    cout << "\nTempo de execucao: " << duration_sec.count() << " segundos" << endl;

    return 0;
}