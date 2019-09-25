#include <iostream>
#include <fstream>
#include <string.h>
#include <bits/stdc++.h>
#include <chrono>

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

void getStartState(vector<vector<int>> &state, vector<vector<int>> &count, float &cost, vector<vector<int>> genes, int k, int v, float **weightMatrix, float CC){
	int len = 0;
	for(int i=0; i<k; i++){
		if(len < genes[i].size()) len = genes[i].size();
	}
	count = vector<vector<int>>(len, vector<int>(v+1, 0));
	for(int i=0; i<k; i++){
		vector<int> vect(len, v);
		for(int j=0; j<len; j++){
			if(j<genes[i].size()) vect[j] = genes[i][j];
			count[j][vect[j]] += 1;	
		}
		state.push_back(vect);
	}
	for(int i=0; i<len; i++){
		for(int l=0; l<v+1; l++){
			for(int m=l; m<v+1; m++){
				cost += count[i][l]*count[i][m]*weightMatrix[l][m];
			}
		}
		cost += count[i][v]*CC;
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

	vector<vector<int>> state;
	vector<vector<int>> count;
	float cost = 0;
	getStartState(state, count, cost, genes, k, vocabSize, weightMatrix, CC);
	//cout <<" COST :" << cost<<"\n";
	vector<int> change(4, -1);
	float update = 0;
	int counter = 0;
	float T = 50;
	float prob;
	float ran;
	float tau = 0.0003;
	int cfactor = cost;
	//cout<<cfactor<<"\n";
	using namespace std::chrono;
	auto start = high_resolution_clock::now(); 
	int limit = 100;
	while(true){
		getNeighbour(state, count, vocabSize, weightMatrix, CC, change);
		update = change[3];
		if(counter == 5){
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start)/5; 
			cfactor = abs(cfactor)*0.85;
			//cout<<cfactor<<"\n";
			//cout << duration.count()<<"\n";
			limit = timeLimit*60*800000/duration.count();
			limit = limit % 700000;
			//cout << limit << endl;
			tau = (float)70/limit;
			//cout<<tau<<endl;
		}
		if(change[3] <= 0){
			updateState(state, change, count, vocabSize);
			cost += update;
		}
		else{
			prob = exp(-update/T);
			ran = (float ) rand() / RAND_MAX;
			//cout<< ran << "  " << prob << "\n";
			if(prob > ran){
				updateState(state, change, count, vocabSize);
				cost += update;
			}
			else{
				//cout<<"---DID NOT UPDATE---"<< "\n";
			}
		}

		
		T = cfactor*exp(-counter*tau);
		if(counter == limit){
			//cout<< ran << "  " << prob << "\n";
			break;
		}
		/*for(int i=0; i<state.size(); i++){
			for(int j = 0; j< state[i].size(); j++){
				cout<<state[i][j]<<" ";
			}
			cout<<"\n";
		}*/
		//cout<<"\n";
		counter++;
	}
	/*for(int i=0; i<state.size(); i++){
		for(int j = 0; j< state[i].size(); j++){
			cout<<vmap[state[i][j]]<<" ";
		}
		cout<<"\n";
	}*/
	//cout<<" Cost : " << cost<<"\n";
	ofstream outFile;
	outFile.open(argv[2]);
	for(int i=0; i<state.size(); i++){
		for(int j = 0; j< state[i].size(); j++){
			outFile<<vmap[state[i][j]];
		}
		outFile<<"\n";
	}
	outFile.close();

}
