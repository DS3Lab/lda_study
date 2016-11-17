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
#define K 40
#define V 12419
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
double Nphi[K][V],Ntheta[M][K];
double nz[K];
int Doc_Count[M];
vector<int> Document[M];
map<string, int> Str2Int;
vector<int> Doc_c[M];
double nzHat[K];
double Gamma[V][K];
double nPhiHat[K][V];
double Kappa = 0.9;

int iteration = 2000;
double S1 = 10,S2 = 50;
double TAO = 1000,TAO2= 105;

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
void  normalize(){

	for (int i = 0;i<M;++i){
		for (int k = 0;k<K;++k){
			Ntheta[i][k] = theta[i][k] + alpha;
		}
	}
	for (int k = 0 ;k<K;++k){
		for (int w = 0;w<V;++w){
			Nphi[k][w] = phi[k][w] + eta;
 		}
	}


	for (int k = 0; k < K; ++k)
	{
		double phisum = 0;
		for (int v = 0; v < V; ++v)
			phisum += Nphi[k][v];

		for (int v = 0; v < V; ++v)
			Nphi[k][v] /= phisum;
	}

	

	for (int i = 0; i < M; ++i)
	{
		double thetasum = 0;

		for (int k = 0; k < K; ++k)
			thetasum += Ntheta[i][k];

		for (int k = 0; k < K; ++k)
			Ntheta[i][k] /= thetasum;
	}

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

				p[k] = Ntheta[m][k] * Nphi[k][wi];
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

void scvb_infer(MiniBatch batch,int iter,int curindex){

	memset(nPhiHat, 0, sizeof(nPhiHat));
	memset(nzHat, 0, sizeof(nzHat));
	memset(Gamma, 0, sizeof(Gamma));
	memset(nz,0,sizeof(nz));
	int doc_iter = 0;
	int burn_Total = 15;
	double rhoTheta ;

	double rhoPhi;


	for (int k = 0;k<K;++k){
		for (int w= 0;w<V;++w){
			nz[k] += phi[k][w];
		}
	}
	double C = GetTotalNum();
	double minibatch_per_corpus =(1.0*C)/batch.BatchSize;

	

	for (int i = 0; i < batch.BatchSize; ++i)
	{
		int index = batch.doc[i].index;
		doc_iter = curindex*batch.BatchSize + iter*M +i+1 ;
		for (int burn = 0; burn < burn_Total; ++burn)
		{

			rhoTheta =  S2/ pow((TAO2+burn ), Kappa);

			double weight_left = 0;
			double GammaK[K];
			memset(GammaK,0,sizeof(GammaK));
			for (int j = 0; j < Document[index].size(); ++j)
			{

				double GammaSum = 0;
				int wi = Document[index][j];
				int wc = Doc_c[index][j];

				for (int k = 0; k < K; ++k)
				{
					Gamma[wi][k] = (eta + phi[k][wi]) * (alpha + theta[index][k]) / (nz[k] + eta * V);
					GammaSum += Gamma[wi][k];
				}
				for (int k = 0;k<K;++k)
					Gamma[wi][k] /= GammaSum;

				weight_left += pow((1 - rhoTheta), wc)   * wc /Doc_Count[index];

				for (int k =0;k<K;++k)
				{
					GammaK[k] += (1 - pow((1 - rhoTheta), wc)) * wc * Gamma[wi][k];
				}
			}		

			for (int k = 0;k<K;++k)
			{
				theta[index][k] = weight_left * theta[index][k] + GammaK[k];

			}
		}

		for (int j = 0; j < Document[index].size(); ++j)
		{
			double GammaSum = 0;
			int wi = Document[index][j];
			int wc = Doc_c[index][j];
			for (int k = 0; k < K; ++k)
			{
				Gamma[wi][k] = (eta + phi[k][wi]) * (alpha + theta[index][k]) / (nz[k] + eta * V);
				GammaSum  += Gamma[wi][k];
			}
			for (int k = 0;k<K;++k)
				Gamma[wi][k] /= GammaSum;

			for(int k = 0;k<K;++k)		{
				double tep = minibatch_per_corpus* Gamma[wi][k] * wc ;

				nPhiHat[k][wi] += tep;
				nzHat[k] += tep ;
			}
		}

	}

	rhoPhi = S1 / pow(TAO +doc_iter, Kappa);

	for (int k = 0; k < K; ++k)
	{
		for (int w = 0; w < V; ++w)
		{
			phi[k][w] = (1 - rhoPhi)*phi[k][w] + rhoPhi * nPhiHat[k][w];
		}
	}

}
vector<MiniBatch> SetBatchs(){
	int BatchSize = 20;
	
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

			scvb_infer(Batchs[i],iter,i);
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
