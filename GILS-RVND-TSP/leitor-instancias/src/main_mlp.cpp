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

struct subsequence{
    double T, C;
    int W, first, last; 
    inline static subsequence concatenate(Data &data, subsequence &sigma_1, subsequence &sigma_2){
        subsequence sigma;
        double temp = data.getDistance(sigma_1.last, sigma_2.first);
        sigma.T = sigma_1.T + temp + sigma_2.T;
        sigma.W = sigma_1.W + sigma_2.W;
        sigma.C = sigma_1.C + sigma_2.W * (sigma_1.T + temp) + sigma_2.C;
        sigma.first = sigma_1.first;
        sigma.last = sigma_2.last;

        return sigma;
    }
};

void updateAllSubseq(Data &data, solution &s, vector<vector<subsequence>> &subseq_matrix){
    int n = s.sequence.size();

    for(int i = 0; i < n; i++){
        int v = s.sequence[i];
        subseq_matrix[i][i].W = (i > 0);
        subseq_matrix[i][i].C = 0;
        subseq_matrix[i][i].T = 0;
        subseq_matrix[i][i].first = s.sequence[i];
        subseq_matrix[i][i].last = s.sequence[i];
    }

    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            subseq_matrix[i][j] = subsequence::concatenate(data, subseq_matrix[i][j-1], subseq_matrix[j][j]);
        }
    }
    for(int i = n - 1; i >= 0; i--){
        for(int j = i - 1; j >= 0; j--){
            subseq_matrix[i][j] = subsequence::concatenate(data, subseq_matrix[i][j+1], subseq_matrix[j][j]);
        }    
    }
}

bool bestImprovementSwap(solution& s, Data& data, vector<vector<subsequence>> &subseq_matrix){
    double bestDelta = 0.0;
    int best_i = 0, best_j = 0;
    int n = s.sequence.size();

    for(int i = 1; i < n - 1; i++){
        for(int j = i + 1; j < n - 1; j++){
            subsequence aux1, aux2, aux3, newSolution;
            double delta;
            
            if(j == i + 1){
                aux1 = subsequence::concatenate(data, subseq_matrix[0][i - 1], subseq_matrix[j][j]);
                aux2 = subsequence::concatenate(data, aux1, subseq_matrix[i][i]);
                newSolution = subsequence::concatenate(data, aux2, subseq_matrix[j + 1][n -1]);
                delta = newSolution.C - s.cost;
            } 
            else{
                aux1 = subsequence::concatenate(data, subseq_matrix[0][i - 1], subseq_matrix[j][j]);
                aux2 = subsequence::concatenate(data, aux1, subseq_matrix[i + 1][j - 1]);
                aux3 = subsequence::concatenate(data, subseq_matrix[i][i], subseq_matrix[j + 1][n - 1]);
                newSolution = subsequence::concatenate(data, aux2, aux3);
                delta = newSolution.C - s.cost;
            }   
            
            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j; 
            }               
        }
    }

     if(bestDelta < 0){
        swap(s.sequence[best_i], s.sequence[best_j]);
        s.cost += bestDelta;

        updateAllSubseq(data, s, subseq_matrix);

        return true;
    }

    return false;
}

bool bestImprovement2Opt(solution& s, Data& data, vector<vector<subsequence>> &subseq_matrix){
    double bestDelta = 0.0;
    int best_i = 0, best_j = 0;
    int n = s.sequence.size();

    for(int i = 1; i < n - 1; i++){
        for(int j = i + 2; j < n - 1; j++){
            subsequence aux1, newSolution;
            double delta;

            aux1 = subsequence::concatenate(data, subseq_matrix[0][i], subseq_matrix[j][i + 1]);
            newSolution = subsequence::concatenate(data, aux1, subseq_matrix[j + 1][n - 1]);
            delta = newSolution.C - s.cost;

            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }               
        }
    }

    if(bestDelta < 0){
        reverse(s.sequence.begin() + best_i + 1, s.sequence.begin() + best_j + 1);
        s.cost += bestDelta;

        updateAllSubseq(data, s, subseq_matrix);
        
        return true;
    }

    return false;
}

bool bestImprovementOrOpt(solution& s, Data& data, int l, vector<vector<subsequence>> &subseq_matrix){
    double bestDelta = 0.0;
    int best_i = 0, best_j = 0;
    int n = s.sequence.size();

    for(int i = 1; i <= n - l - 1; i++){
        for(int j = 1; j <= n - l - 1; j++){
            subsequence aux1, aux2, newSolution;
            double delta;

            // Verifica sobreposição
            if(j >= i && j <= i + l - 1) continue;
            
            // Pula movimentos sem efeito (adjacentes)
            if(j == i - 1 || j == i + l) continue;
            
            if(j > i){
                aux1 = subsequence::concatenate(data, subseq_matrix[0][i - 1], subseq_matrix[j][j + l -1]);
                aux2 = subsequence::concatenate(data, aux1, subseq_matrix[i][j - 1]);
                newSolution = subsequence::concatenate(data, aux2, subseq_matrix[j + l][n - 1]);
                delta = newSolution.C - s.cost;
            }
            else{
                aux1 = subsequence::concatenate(data, subseq_matrix[0][j - 1], subseq_matrix[j + l][i]);
                aux2 = subsequence::concatenate(data, aux1, subseq_matrix[j][j + l -1]);
                newSolution = subsequence::concatenate(data, aux2, subseq_matrix[i + 1][n - 1]);
                delta = newSolution.C - s.cost;
            }    

            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }               
        }
    }

    if(bestDelta < 0){
        vector<int> bloco;
        for(int k = 0; k < l; k++){
            bloco.push_back(s.sequence[best_i + k]);
        }

        s.sequence.erase(s.sequence.begin() + best_i, s.sequence.begin() + best_i + l);
        
        int posInsercao;
        if(best_j < best_i){
            posInsercao = best_j;      // Inserir antes
        } 
        else {
            posInsercao = best_j - l;  // Inserir depois
        }
        
        s.sequence.insert(s.sequence.begin() + posInsercao, bloco.begin(), bloco.end());
        s.cost += bestDelta;

        updateAllSubseq(data, s, subseq_matrix);

        return true;
    }

    return false;
}

void buscaLocal(solution& s, Data& data, vector<vector<subsequence>> &subseq_matrix){
    vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;

    while(!NL.empty()){
        int v = rand() % NL.size();

        switch(NL[v]){
            case 1:
                improved = bestImprovementSwap(s, data, subseq_matrix); 
                break;
            case 2:
                improved = bestImprovement2Opt(s, data, subseq_matrix); 
                break;
            case 3:
                improved = bestImprovementOrOpt(s, data, 1, subseq_matrix); 
                break;
            case 4:
                improved = bestImprovementOrOpt(s, data, 2, subseq_matrix);
                break;
            case 5:
                improved = bestImprovementOrOpt(s, data, 3, subseq_matrix);
                break;           
        }

        if(improved){
            NL = {1, 2, 3, 4, 5};
        } else {
            NL.erase(NL.begin() + v);
        }
    }
}

solution perturbacao(solution best, Data& data){
    solution copia = best;
    int n = copia.sequence.size() - 1;
    
    int maxBlock = max(2, (int)ceil(n / 10.0));
    int tamanho1 = 2 + rand() % (maxBlock - 1);
    int tamanho2 = 2 + rand() % (maxBlock - 1);
    
    int maxPos1 = n - tamanho1 - tamanho2 - 3;
    int pos1 = 1 + rand() % maxPos1;
    
    int minPos2 = pos1 + tamanho1 + 1;
    int maxPos2 = n - tamanho2 - 1;
    int pos2 = minPos2 + rand() % (maxPos2 - minPos2 + 1);
    
    vector<int> bloco1(copia.sequence.begin() + pos1, copia.sequence.begin() + pos1 + tamanho1);
    vector<int> bloco2(copia.sequence.begin() + pos2, copia.sequence.begin() + pos2 + tamanho2);
    
    vector<int> novaSequencia;
    novaSequencia.insert(novaSequencia.end(), copia.sequence.begin(), copia.sequence.begin() + pos1);
    novaSequencia.insert(novaSequencia.end(), bloco2.begin(), bloco2.end());
    novaSequencia.insert(novaSequencia.end(), copia.sequence.begin() + pos1 + tamanho1, copia.sequence.begin() + pos2);
    novaSequencia.insert(novaSequencia.end(), bloco1.begin(), bloco1.end());
    novaSequencia.insert(novaSequencia.end(), copia.sequence.begin() + pos2 + tamanho2, copia.sequence.end());
    
    copia.sequence = novaSequencia;

    return copia;
}

solution ILS(int maxIter, int maxIterIls, Data& data){
    solution bestOfAll;
    bestOfAll.cost = 1e9;

    for(int i = 0; i < maxIter; i++){
        int iterIls = 0;

        double alpha = (double)rand() / RAND_MAX;
        solution s = construcao(data, alpha);
            
        vector<vector<subsequence>> subseq_matrix(s.sequence.size(), vector<subsequence>(s.sequence.size()));

        updateAllSubseq(data, s, subseq_matrix);
        s.cost = subseq_matrix[0][s.sequence.size() - 1].C;

        solution best = s;

        while(iterIls < maxIterIls){
            buscaLocal(s, data, subseq_matrix);
            
            if(s.cost < best.cost){
                best = s;
                iterIls = 0;
            }

            s = perturbacao(best, data);
            updateAllSubseq(data, s, subseq_matrix);
            s.cost = subseq_matrix[0][s.sequence.size() - 1].C;
            
            iterIls++;
        }

        if (best.cost < bestOfAll.cost){
            bestOfAll = best;
        }
    }

    return bestOfAll;
}

int main(int argc, char** argv){
    auto start = chrono::high_resolution_clock::now();

    srand(time(NULL));

    auto data = Data(argc, argv[1]);
    data.read();
    size_t n = data.getDimension();
    
    int maxIter = 50;
    int maxIterILS;
    if (n >= 150){
        maxIterILS = n / 2;
    } else{
        maxIterILS = n;
    }

    cout << "Dimension: " << n << endl;

    cout << "\nExecutando ILS para " << n << " cidades..." << endl;

    solution bestSolution = ILS(maxIter, maxIterILS, data);

    cout << "Melhor solucao encontrada: 0 -> ";
    for (int i = 0; i < (int)n; i++){
        cout << bestSolution.sequence[i] << " -> ";
    }

    cout << bestSolution.sequence[0] << " -> 0" << endl;

    cout << "Custo da melhor solucao: " << bestSolution.cost << endl;

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration_sec = end - start;
    cout << "\nTempo de execucao: " << duration_sec.count() << " segundos" << endl;

    return 0;
}