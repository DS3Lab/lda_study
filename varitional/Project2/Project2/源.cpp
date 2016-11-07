#include<iostream>
#include<stdio.h>
#include<string>
#include<fstream>
#include<sstream>
#include<map>
#include<vector>
#include<cmath>
#include<algorithm>
#define M 1000 
#define K 40
#define V 11980
#define MaxRand 32761
#define MaxLine 200000
using namespace std;



int nmk[M][K], nkt[K][V], nktS[K];
vector<int> Z[M];


double alpha = 0.1, beta = 0.01;
double phi[K][V], theta[M][K];
vector<int> Document[M];
map<string, int> Str2Int;
vector<int> Doc_c[M];

vector<vector<double>> Gamma[M];
double m_mk[M][K], v_mk[M][K], m_kt[K][V], v_kt[K][V];
double m_k[K], v_k[K];
int iteration = 500;


void count_update(int doc, int j, int word, double scale){
	for (int k = 0; k < K; ++k){
		double mc = scale * Gamma[doc][j][k];
		double vc = scale * (1 - Gamma[doc][j][k]) * (Gamma[doc][j][k]);
		m_mk[doc][k] += mc;
		v_mk[doc][k] += vc;
		m_kt[k][word] += mc;
		v_kt[k][word] += vc;
		m_k[k] += mc;
		v_k[k] += vc;
	}
	return;
}


void GetCount(){

	for (int i = 0; i < M; ++i)
		Doc_c[i].resize(Document[i].size());
	for (int i = 0; i < M; ++i)
		for (int j = 0; j < Doc_c[i].size(); ++j)
			Doc_c[i][j] = 0;


	for (int i = 0; i < M; ++i)
	{
		int index = -1;
		for (int j = 0; j < Document[i].size(); ++j)
		{
			if (j - 1 >= 0 && Document[i][j] == Document[i][j - 1])
			{
				Doc_c[i][index]++;

			}
			else if (j - 1 < 0 || (j - 1 >= 0 && Document[i][j] != Document[i][j - 1]))
			{
				index++;
				Doc_c[i][index]++;
			}
		}
		Doc_c[i].resize(index + 1);
	}

}

void GetDocment(string path){
	ifstream in(path, ios::in);//以输入方式打开文件 
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
		for (auto it : vec){
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
		
		Document[i].erase(p,Document[i].end());
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
		for (int j = 0; j < K; ++j)
			theta[i][j] = 0;
	for (int i = 0; i < K; ++i)
		for (int j = 0; j < V; ++j)
			phi[i][j] = 0;
	memset(m_mk, 0, sizeof(m_mk));
	memset(v_mk, 0, sizeof(v_mk));
	memset(m_kt, 0, sizeof(m_kt));
	memset(v_kt, 0, sizeof(v_kt));
	memset(m_k, 0, sizeof(m_k));
	memset(v_k, 0, sizeof(v_k));

	for (int i = 0; i < M; ++i)
	{
		Gamma[i].resize(Document[i].size());
		for (int j = 0; j < Document[i].size(); ++j)
			Gamma[i][j].resize(K);
	}

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < Document[i].size(); ++j)
		{
			double teptotal = 0;
			for (int k = 0; k < K; ++k)
			{
				Gamma[i][j][k] = getrand();
				teptotal += Gamma[i][j][k];
			}
			
			//cout << teptotal << endl;
			/*
			for (int k = 0; k < K; ++k)
			{
				Gamma[i][j][k] = Gamma[i][j][k] / (teptotal);
				
				double temp = Doc_c[i][j] * Gamma[i][Document[i][j]][k];
				theta[i][k] += temp;
				phi[Document[i][j]][k] += temp;
				
			}*/
		}
	}
	


	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < Document[i].size(); ++j)
		{
			count_update(i, j, Document[i][j], Doc_c[i][j]);
		}
	}

	cout << "init end" << endl;
	return;

}
double llhw(){
	double lgamma_alpha = lgamma(alpha);
	double lgamma_beta = lgamma(beta);



	double res = 0;
	res += K* lgamma(beta*V);
	for (int k = 0; k < K; ++k)
		res -= lgamma(beta*V + m_k[k]);



	for (int t = 0; t < V; ++t){
		for (int k = 0; k < K; ++k){
			if (m_kt[k][t] > 0)
				res += lgamma(beta + m_kt[k][t]) - lgamma_beta;
		}
	}

	for (int m = 0; m < M; ++m)
	{
		res += (lgamma(alpha*K) - lgamma(alpha*K + Document[m].size()));
		for (int k = 0; k < K; ++k){
			if (m_mk[m][k] > 0)
				res += lgamma(alpha + m_mk[m][k]) - lgamma_alpha;
		}

	}
	return res;
}


double cvb_infer(int iter){
	int D = M;
	double maxdelta_Gamma = 0;
	double max_delta = 0;

	for (int i = 0; i < D; ++i){

		for (int j = 0; j < Document[i].size(); ++j)
		{

			int wi = Document[i][j];
			int ci = Doc_c[i][j];

			double normal = 0;
			count_update(i, j, wi, -ci);

			double old_Gamma[K];
			for (int k = 0; k < K; ++k)
			{

				old_Gamma[k] = Gamma[i][j][k];
				double newGamma = (alpha + m_mk[i][k]) *
					((beta + m_kt[k][wi]) / (beta*V + m_k[k]))*
					
					exp( -1 * (v_mk[i][k] / (2 * pow(alpha + m_mk[i][k], 2))) -
					(v_kt[k][wi] / (2 * pow(beta + m_kt[k][wi], 2))) +
					(v_k[k] / (2 * pow(beta*V + m_k[k], 2))));
				Gamma[i][j][k] = newGamma;
				normal += newGamma;
			}

			double cur_delta = 0;
			//cout << normal << endl;
			for (int k = 0; k < K; ++k){
				Gamma[i][j][k] = Gamma[i][j][k] / (normal);

				cur_delta += abs((Gamma[i][j][k] - old_Gamma[k]));
			}
			if (cur_delta > max_delta)
				max_delta = cur_delta;

			count_update(i, j, wi, ci);

		}
	}
	return max_delta;

}

//get phi from Gamma
void est_theta(){
	double m_msum[M];
	for (int i = 0; i < M; ++i)
	{
		m_msum[i] = 0;
		for (int j = 0; j < K; ++j)
		{
			m_msum[i] += m_mk[i][j];
		}

	}

	for (int i = 0; i < M; ++i)
	{
		for (int k = 0; k < K; ++k){
			theta[i][k] = (alpha + m_mk[i][k]) / (alpha*K + m_msum[i]);
		}
	}


}


// get theta from theta

void est_phi(){
	double m_vsum[K];

	for (int i = 0; i < K; ++i)
	{
		m_vsum[i] = 0;
		for (int v = 0; v < V; ++v)
		{
			m_vsum[i] += m_kt[i][v];
		}
	}
	for (int i = 0; i < K; ++i)
	{
		for (int v = 0; v < V; ++v)
		{

				phi[i][v] = (beta + m_kt[i][v]) / (beta*V + m_vsum[i]);

		}
	}


}
void PrintTheta(string path){
	fstream out(path, ios::out);
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < K; ++j)
			out << theta[i][j] << " ";
		out << endl;
	}
}

int  GetTotalNum(){
	int res = 0;
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < Doc_c[i].size();++j)
		res += Doc_c[i][j];
	}
	return res;

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
			result += Doc_c[m][j]*log(tepsum);
		}

	}
	int TotalNum = GetTotalNum();
	result = (-result / TotalNum);

	return result;
}

void Train(){

	for (int iter = 0; iter < iteration; ++iter){
		cout << "iter  " << iter << endl;
		double infer = llhw();
		double tepres = cvb_infer(iter);

		cout << "max_change" << tepres << endl;
		cout << "llhw" << infer << endl;
		est_phi();
		est_theta();
		double perplexity = Compute_Perplexity();
		cout << "perplexity  " << perplexity << endl;
	}
	return;

}
void PrintResult(string path){
	ofstream out(path);
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
	int t = 12;
	srand(t);
	GetDocment("test.txt");
	GetCount();
	cout << "Get count end!!" << endl;
	RemoveDocment();
	Init();

	Train();

	est_phi();
	est_theta();
	PrintTheta("theta.txt");
	PrintResult("result2.txt");
	return 0;

}