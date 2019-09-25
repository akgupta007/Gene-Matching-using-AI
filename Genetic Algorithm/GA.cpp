#include <iostream>
#include <fstream>
#include <string.h>
#include <bits/stdc++.h>

using namespace std;

// adjust your algorithm with time


void parse(char* name, float &timeLimit, int &vocabSize, map<char, int> &vocab, int &k, vector<vector<int>> &genes, float &CC, float** &weightMatrix){
	ifstream inFile;
	inFile.open(name);
	if(!inFile){
		cout<<"Error opening File";
		exit(0);
	}
	string str;

	getline(inFile, str);
	timeLimit = stof(str);
	
	getline(inFile, str);
	vocabSize = stoi(str);
	getline(inFile, str);
	
	k =0;
	for(int i=0; i<str.length(); i++){
		if(str[i] > 64 && str[i] <91){
			vocab.insert({str[i], k++});
		}
	}
	getline(inFile,str);
	k = stoi(str);
	int l;
	for(int i =0 ; i<k; i++){
		getline(inFile, str);
		l = str.length();
		if(str.compare(l-1, 1, " ")){
			l--;
		}
		vector<int> vect(l, 0);
		genes.push_back(vect);
		for(int j=0; j< l; j++){
			genes[i][j] = vocab[str[j]];
		}
		str.clear();
	}
	getline(inFile, str);
	CC = stof(str);
	
	int i =0; 
	int j =0;
	weightMatrix = new float*[vocabSize+1];
	for(int i =0; i<vocabSize+1; i++){
		weightMatrix[i] = new float[vocabSize+1];
		getline(inFile, str);
		stringstream s(str);
		for(int j= 0; j<vocabSize+1; j++){
			s >> weightMatrix[i][j];
		}
	}
	inFile.close();
}

void getNeighbour(vector<vector<int>> &state, vector<vector<int>> count, int v, float** weightMatrix, float CC, vector<int> &change){
	
	vector<vector<int>> allNeighbours;

	// backward pass
	float rShiftCost = 0;
	int shiftState = 0;
	float temp_cost;
	int cmin = 2147483647;
	int rvalue = -1;
	for(int i = 0; i<state.size(); i++){
		rShiftCost = 0;
		shiftState = 0;
		temp_cost = 0;
		rvalue = -1;
		for(int j = state[i].size()-1; j>=0; j--){

			if(shiftState == 0 && state[i][j] != v){
				shiftState = 1; 
				rvalue = j+1;
			}
			if(shiftState == 1 && state[i][j] != v){

				if(j+1 >= state[i].size()){
					rShiftCost = (state.size()-1)*weightMatrix[state[i][j]][v] + state.size()*CC;
				}
				else{
					count[j+1][state[i][j+1]] -= 1;
					for(int l = 0; l<v+1;l++){
						rShiftCost = rShiftCost - count[j+1][l]*weightMatrix[state[i][j+1]][l] + count[j+1][l]*weightMatrix[state[i][j]][l];
					}
					count[j+1][state[i][j+1]] += 1;
				}
				//cout<<rShiftCost <<" ";


				temp_cost = 0;
				count[j][state[i][j]] -= 1;
				for(int l = 0; l<v+1;l++){
					temp_cost = temp_cost - count[j][l]*weightMatrix[state[i][j]][l] + count[j][l]*weightMatrix[v][l];
				}
				if(count[j][v] == state.size()-1){
					temp_cost -= state.size()*CC;
				}
				count[j][state[i][j]] += 1;
				cmin = (int ) rShiftCost + temp_cost;
				allNeighbours.push_back(vector<int>{i, j, rvalue, cmin});
				
			}
			if(state[i][j] == v){
				shiftState = 0;
				rShiftCost = 0;
				rvalue = -1;
			}
		}
		//cout<<"\n";
	}

	// forward pass
	float lShiftCost = 0;
	int lvalue = -1;
	for(int i = 0; i<state.size(); i++){
		shiftState = -1;
		lShiftCost = 0;
		lvalue = -1;
		temp_cost = 0;
		for(int j = 0; j < state[i].size(); j++){
			
			if(shiftState == 0 && state[i][j] != v){
				shiftState = 1; 
				lvalue = j-1;
			}

			if(shiftState == 1 && state[i][j] != v){

				count[j-1][state[i][j-1]] -= 1;
				for(int l = 0; l<v+1;l++){
					lShiftCost = lShiftCost - count[j-1][l]*weightMatrix[state[i][j-1]][l] + count[j-1][l]*weightMatrix[state[i][j]][l];
				}
				count[j-1][state[i][j-1]] += 1;
			

				temp_cost = 0;
				count[j][state[i][j]] -= 1;
				for(int l = 0; l<v+1;l++){
					temp_cost = temp_cost - count[j][l]*weightMatrix[state[i][j]][l] + count[j][l]*weightMatrix[v][l];
				}
				if(count[j][v] == state.size()-1){
					temp_cost -= state.size()*CC;
				}
				count[j][state[i][j]] += 1;
				cmin =  lShiftCost + temp_cost;

				allNeighbours.push_back(vector<int>{i, lvalue, j, cmin});
			}
			if(state[i][j] == v){
				shiftState = 0;
				lShiftCost = 0;
				lvalue = -1;
			}
		}
	}




	///////// remove dashes in the  :: state --> DONE

	int n = rand() % allNeighbours.size();
	change = allNeighbours[n];
}

void getStartState(vector<vector<int>>* &state, vector<vector<int>>* &count, int popSize, vector<int> &cost, vector<vector<int>> genes, int k, int v, float **weightMatrix, float CC){
	int len = 0;
	for(int i=0; i<k; i++){
		if(len < genes[i].size()) len = genes[i].size();
	}
	count[0] = vector<vector<int>>(len, vector<int>(v+1, 0));
	for(int i=0; i<k; i++){
		vector<int> vect(len, v);
		for(int j=0; j<len; j++){
			if(j<genes[i].size()) vect[j] = genes[i][j];
			count[0][j][vect[j]] += 1;	
		}
		state[0].push_back(vect);
	}
	int lg;
	for(int p=1; p<popSize;p++){
		state[p] = state[0];
		count[p] = vector<vector<int>>(len, vector<int>(v+1, 0));
		for(int i=0 ;i<k; i++){
			lg=0;
			std::random_shuffle ( state[p][i].begin(), state[p][i].end() );
			for(int j=0; j<state[p][i].size();j++){
				if(state[p][i][j] != v){
					state[p][i][j] = genes[i][lg++];
				}
				count[p][j][state[p][i][j]] += 1;
			}
		}
	}
	for(int p=0; p<popSize;p++){
		for(int i=0; i<len; i++){
			for(int l=0; l<v+1; l++){
				for(int m=l; m<v+1; m++){
					cost[p] += count[p][i][l]*count[p][i][m]*weightMatrix[l][m];
				}
			}
			cost[p] += count[p][i][v]*CC;
		}
	}
}

void updateState(vector<vector<int>> &state, vector<int> &change,vector<vector<int>> &count, int v){
	int index = change[0];
	int lvalue = change[1];
	int rvalue = change[2];
	if(state[index][lvalue] == v){
		for(int i = lvalue; i<rvalue; i++){
			count[i][state[index][i]] -= 1;
			count[i][state[index][i+1]] += 1;
			state[index][i] = state[index][i+1];
		}
		count[rvalue][state[index][rvalue]] -= 1;
		count[rvalue][v] += 1;
		state[index][rvalue] = v;

		if(count[rvalue][v] == state.size()){
			for(int i=0; i<state.size(); i++){
				state[i].erase(state[i].begin() + rvalue);
			}
			count.erase(count.begin() + rvalue);
		}

	}
	else{
		if(rvalue >= state[0].size()){
			for(int i=0; i<state.size(); i++){
				state[i].push_back(v);
			}
			count.push_back(vector<int>(v+1,0));
			count[rvalue][v] = state.size();
		}
		for(int i=rvalue; i>lvalue; i--){
			count[i][state[index][i]] -= 1;
			count[i][state[index][i-1]] += 1;
			state[index][i] = state[index][i-1];
		}
		count[lvalue][state[index][lvalue]] -= 1;
		count[lvalue][v] += 1;
		state[index][lvalue] = v;
		if(count[lvalue][v] == state.size()){
			for(int i=0; i<state.size(); i++){
				state[i].erase(state[i].begin() + lvalue);
			}
			count.erase(count.begin() + lvalue);
		}
	}
}



void crossOver(vector<vector<int>>* &state, vector<int> &cost, vector<vector<int>> *count, int n1[], int n2[], int popSize, int v, float CC, float** weightMatrix){
	vector<vector<int>> *newState = new vector<vector<int>>[popSize];
	vector<vector<int>> *newCount = new vector<vector<int>>[popSize];
	for(int p=0; p<popSize; p++){
		int r = rand() % state[n1[p]][0].size()/3+ state[n1[p]][0].size()/3;
		int cnt[state[0].size()] = {0};
		for(int i=0; i<state[n1[p]].size(); i++){
			for(int j=0; j<r; j++){
				if(state[n1[p]][i][j] != v){
					cnt[i]++;
				}
			}
		}
		int max = 0;

		for(int i=0; i<state[n2[p]].size(); i++){
			int c = 0;
			for(int j=0; j<state[n2[p]][0].size(); j++){
				if(cnt[i] == c && state[n2[p]][i][j] != v){
					int ns = state[n2[p]][0].size() - j;
					if(max < ns){
						max = ns;
					}
					break;
				}
				if(state[n2[p]][i][j] != v){
					c++;
				}
			}
			
		}

		int nsize = r+max;
		
		newCount[p] = vector<vector<int>>(nsize, vector<int>(v+1, 0));
		newState[p] = vector<vector<int>>(state[n1[p]].size(), vector<int>(nsize, v));
		for(int i=0; i<state[n1[p]].size(); i++){
			for(int j=0; j<r; j++){
				newState[p][i][j] = state[p][i][j];
				newCount[p][j][newState[p][i][j]] += 1;
			}
		}
		int temp = state[n2[p]][0].size()-1;
		for(int i=0; i<state[n2[p]].size(); i++){
			for(int j=0; j<max; j++){
				if(cnt[i] != 0){
					if(state[p][i][temp-j] != v){
						cnt[i]--;
						newState[p][i][nsize - j -1 ] = state[p][i][temp-j];
					}
				}
				
				newCount[p][nsize - j -1 ][newState[p][i][nsize - j -1 ]] += 1;
			}
		}
	}

	/////remove vacant line

	state = newState;
	count = newCount;
	for(int p=0; p<popSize;p++){
		cost[p] = 0;
		for(int i=0; i<state[p][0].size(); i++){
			for(int l=0; l<v+1; l++){
				for(int m=l; m<v+1; m++){
					cost[p] += count[p][i][l]*count[p][i][m]*weightMatrix[l][m];
				}
			}
			cost[p] += count[p][i][v]*CC;
		}
	}
	
}

template <class T1, class T2>
map<T2, T1> swapPairs(map<T1, T2> m) {
    map<T2, T1> m1;

    for (auto&& item : m) {
        m1.emplace(item.second, item.first);
    }

    return m1;
};


int main(int argc, char** argv){
	float timeLimit, CC;
	int k, vocabSize;
	map<char, int> vocab;
	vector<vector<int>> genes;
	float** weightMatrix;
	parse(argv[1], timeLimit, vocabSize, vocab, k, genes, CC, weightMatrix);
	vocab.insert({'-', vocabSize});
	map<int,char> vmap = swapPairs(vocab);

	int popSize = 10;
	vector<vector<int>> *state = new vector<vector<int>>[popSize];
	vector<vector<int>> *count = new vector<vector<int>>[popSize];
	vector<int> cost(popSize, 0);
	getStartState(state, count, popSize, cost, genes, k, vocabSize, weightMatrix, CC);
	for(int p= 0; p<popSize; p++){
		cout <<" COST :" << cost[p]<<"\n";
		for(int i=0; i<state[p].size(); i++){
			for(int j = 0; j< state[p][i].size(); j++){
				cout<<vmap[state[p][i][j]]<<" ";
			}
			cout<<"\n";
		}
		cout<<"\n";

	}
	int n1[popSize];
 	int n2[popSize];
 	int counter = 0;
 	while(counter<4){
 		cout<<counter << " : \n";
 		default_random_engine generator;
	 	discrete_distribution<int> distribution(cost.begin(), cost.end());
	 	for(int i=0; i<popSize; i++){
	 		n1[i] = distribution(generator);
	 		n2[i] = distribution(generator);
	 		while(n1[i] == n2[i]){
	 			n2[i] = distribution(generator);
	 		}
	 	}
	 	crossOver(state, cost, count, n1, n2, popSize, vocabSize, CC, weightMatrix);
	 	for(int p= 0; p<popSize; p++){
			cout <<" COST :" << cost[p]<<"\n";
			for(int i=0; i<state[p].size(); i++){
				for(int j = 0; j< state[p][i].size(); j++){
					cout<<vmap[state[p][i][j]]<<" ";
				}
				cout<<"\n";
			}
			cout<<"\n";

		}
		cout<<" \n\n";
		counter++;
 	}

}
