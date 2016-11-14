#include<iostream>
#include<stdio.h>
#include<string>
#include<fstream>
#include<sstream>
#include<map>
#include<vector>
#include<cmath>
#include<algorithm>
#include<string.h>
#include<memory.h>
#define M 1500 
#define K 80
#define V 12387
#define MaxRand 32761
#define MaxLine 200000
using namespace std;



class Doc{
public :
	int index;
	int TokenNum;
	vector<int> words;

};
class MiniBatch{
public:
	vector<Doc> doc;
	int TokenNum;
	int BatchSize;

};
double alpha = 0.1, eta = 0.01;
double phi[K][V], theta[M][K];
double nz[K];
int Doc_Count[M];
int TotalNum = 0;
vector<int> Document[M];
map<string, int> Str2Int;
vector<int> Doc_c[M];
double nzHat[K];
double Gamma[V][K];
double nPhiHat[K][V];
double s = 1;
double Kappa = 0.9;

int iteration = 2000;
int TotalTokens = 0;


void GetCount(){

	for (int i = 0; i < M; ++i)
		Doc_c[i].resize(Document[i].size());

	for (int i = 0; i < M; ++i)
		for (int j = 0; j < Doc_c[i].size(); ++j)
			Doc_c[i][j] = 0;

	memset(Doc_Count, 0, sizeof(Doc_Count));

	for (int i = 0; i < M; ++i)
	{
		int index = -1;
		for (int j = 0; j < Document[i].size(); ++j)
		{
			if (j - 1 >= 0 && Document[i][j] == Document[i][j - 1])
			{
				Doc_c[i][index]++;
			
			}
			else if (j-1 < 0 || (j-1>=0 && Document[i][j] != Document[i][j-1]))
			{
				index++;
				Doc_c[i][index]++;
			}
		}
		Doc_c[i].resize(index+1);
	}

	for (int i = 0; i < M; ++i)
		for (int j = 0; j < Document[i].size(); ++j)
			Doc_Count[i] = Document[i].size();

	for (int i = 0; i < M; ++i)
		TotalNum += Doc_Count[i];


}

void GetDocment(string path){
	ifstream in(path.c_str()); 
	int item;
	char line[MaxLine];
	int docidx = 0;
	int ItemInt = 0;
	while (!in.eof())
	{
		memset(line, 0, sizeof(line));

		in.getline(line, sizeof(line));
		int len = strlen(line);
		string item = "";
		istringstream buffin(line);

		vector<string> vec;

		while (!buffin.eof())
		{
			buffin >> item;
			vec.push_back(item);
		}
		Document[docidx].resize(vec.size());
		int idx = 0;
		for (int d = 0 ; d< vec.size();++d){
			string it = vec[d];
			if (!Str2Int[it])
				Str2Int[it] = ItemInt++;
			Document[docidx][idx++] = Str2Int[it];
		}
		docidx++;
		cout << docidx << endl;
	}
	cout << ItemInt << endl;
}
void RemoveDocment(){
	for (int i = 0; i < M; ++i)
	{
		sort(Document[i].begin(), Document[i].end());
		vector<int>::iterator p = unique(Document[i].begin(), Document[i].end());
		//cout << "original" << i << Document[i].size()<< endl;

		Document[i].erase(p, Document[i].end());
		//cout << "new" <<i<< Document[i].size() << endl;
	}
	cout << "Remove Document end !!" << endl;
}
double getrand(){
	double res = rand();
	return res / RAND_MAX;
}
void Init(){

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < K; ++j)
		{
			theta[i][j] = getrand();
			nz[j] += theta[i][j];
		}
	}
	for (int j = 0; j < V; ++j)
	{
		for (int i = 0; i < K; ++i)
		{
			phi[i][j] = getrand();
		}
	}


	cout << "init end" << endl;
	return;

}

MiniBatch GetMiniBatch(int start, int end){
	MiniBatch result;
	result.TokenNum = 0;
	
	for (int i = start; i < end; ++i){
		Doc tep;
		tep.TokenNum = 0;
		tep.index = i;
		tep.TokenNum = Doc_Count[i];

		for (int j = 0; j < Document[i].size(); ++j){
			tep.words.push_back(Document[i][j]);
			
		}
		result.doc.push_back(tep);
		result.TokenNum += tep.TokenNum;
	}
	return result;
}
int GetTotalNum(){
	int result = 0;
	for (int i = 0; i < M; ++i)
		result += Doc_Count[i];
	return result;
}
void PrintTheta(string path){
	fstream out(path.c_str());
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < K; ++j)
			out << theta[i][j] << " ";
		out << endl;
	}
}

void scvb_infer(MiniBatch batch){

	int rhoPhi_t = 1;
	int rhoTheta_t = 1;
	TotalTokens = GetTotalNum();
	memset(nPhiHat, 0, sizeof(nPhiHat));
	memset(nzHat, 0, sizeof(nzHat));
	memset(Gamma, 0, sizeof(Gamma));


	int burn_Total = 1;
	double rhoTheta = s / pow( (10 + rhoTheta_t), Kappa);

	double rhoPhi = s / pow((10  + rhoPhi_t), Kappa);
	for (int i = 0; i < batch.BatchSize; ++i)
	{
		int index = batch.doc[i].index;

		for (int burn = 0; burn < burn_Total; ++burn)
		{
			rhoTheta = s / pow((10 + rhoTheta_t), Kappa);
			rhoTheta_t++;

			for (int j = 0; j < Document[index].size(); ++j)
			{
				double GammaSum = 0;
				int wi = Document[index][j];
				for (int k = 0; k < K; ++k)
				{
					Gamma[wi][k] = (eta + phi[k][wi]) * (alpha + theta[index][k]) / (nz[k] + eta * batch.TokenNum);

					double clumpconst = pow((1 - rhoTheta), Doc_c[index][j]);
					double partone = clumpconst * theta[index][k];
					double parttwo = (1 - clumpconst) *Doc_Count[index] * Gamma[wi][k];
					theta[index][k] = partone + parttwo;
				}
			}
		}

		rhoTheta = s / pow((10 + rhoTheta_t), Kappa);
		rhoTheta_t++;
		for (int j = 0; j < Document[index].size(); ++j)
		{
			double GammaSum = 0;
			int wi = Document[index][j];

			for (int k = 0; k < K; ++k)
			{
				Gamma[wi][k] = (eta + phi[k][wi]) * (alpha + theta[index][k]) / (nz[k] + eta * batch.TokenNum);
				
				double clumpconst = pow((1 - rhoTheta), Doc_c[index][j]);
				double partone = clumpconst * theta[index][k];
				double parttwo = (1 - clumpconst) *Doc_Count[index] * Gamma[wi][k]  ;
				theta[index][k] = partone + parttwo;
				nPhiHat[k][wi] = nPhiHat[k][wi] + (TotalNum* Gamma[wi][k]/batch.TokenNum  );
				
				nzHat[k] = nzHat[k] + (TotalNum*Gamma[wi][k]/batch.TokenNum );
			}
		}
	}

	rhoPhi = s / pow(100 + rhoPhi_t, Kappa);
	rhoPhi_t++;
	for (int k = 0; k < K; ++k)
	{
		for (int w = 0; w < V; ++w)
		{
			phi[k][w] = (1 - rhoPhi)*phi[k][w] + rhoPhi * nPhiHat[k][w];
		}
		nz[k] = (1 - rhoPhi) *nz[k] + rhoPhi*nzHat[k];
	}

}
void  normalize(){

	for (int k = 0; k < K; ++k)
	{
		double phisum = 0;

		for (int v = 0; v < V; ++v)
		{
			phisum += phi[k][v];

		}
		for (int v = 0; v < V; ++v)
			phi[k][v] /= phisum;
	}

	

	for (int i = 0; i < M; ++i)
	{
		double thetasum = 0;

		for (int k = 0; k < K; ++k)
		{
			thetasum += theta[i][k];
		
		}
		for (int k = 0; k < K; ++k)
			theta[i][k] /= thetasum;
	}

	cout << theta[0][0] << endl;
}
double  Compute_Perplexity(){
	double  result = 0;
	for (int m = 0; m < M; ++m)
	{
		for (int j = 0; j < Document[m].size(); ++j)
		{
			int wi = Document[m][j];
			double p[K];
			for (int i = 0; i < K; ++i)
				p[i] = 0;

			for (int k = 0; k < K; ++k){

				p[k] += theta[m][k] * phi[k][wi];
			}
			double tepsum = 0;

			for (int k = 0; k < K; ++k)
				tepsum += p[k];
			result += Doc_c[m][j] * log(tepsum);
		}

	}
	int TotalNum = GetTotalNum();
	result = (-result / TotalNum);

	return exp(result);
}
vector<MiniBatch> SetBatchs(){
	int BatchSize = 100;
	
	vector<MiniBatch> Batchs;
	int index = 0;
	while (index < M){
		int start = index;
		int end = min(index + BatchSize,M);
		MiniBatch curres = GetMiniBatch(start, end);
		curres.BatchSize = end - start;
		Batchs.push_back(curres);
		index = end;
	}
	return Batchs;

}
void Train(){
	vector<MiniBatch> Batchs = SetBatchs();

	for (int iter = 0; iter < iteration; ++iter){
		cout << "iter  " << iter << endl;
		for (int i = 0; i < Batchs.size(); ++i){

			scvb_infer(Batchs[i]);
		}
		normalize();
		double perplexity = Compute_Perplexity();

		cout << "perplexity  " << perplexity << endl;
	}

	return;

}
void outputtest(){
	


}
void PrintResult(string path){
	ofstream out(path.c_str());
	for (int i = 0; i < M; ++i)
	{
		int Maxindex = 0; double MaxValue = theta[i][0];

		for (int k = 0; k < K; ++k)
		{
			if (theta[i][k] > MaxValue){
				MaxValue = theta[i][k];
				Maxindex = k;
			}
		
		}
		out << Maxindex << endl;
	
	}
	out.close();
}
int main(){

	GetDocment("nips.train.txt");
	GetCount();
	cout << "Get count end!!" << endl;

	RemoveDocment();
	Init();

	Train();


	PrintTheta("theta.txt");
	PrintResult("result.txt");
	return 0;

}
