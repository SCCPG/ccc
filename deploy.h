#ifndef __DEPLOY_H__
#define __DEPLOY_H__

#include "lib_io.h"
////////////////////////////////
#include<queue>
#include<iostream>
#include<algorithm>
#include<vector>
#include<bitset>
#include<ctime>
#include<cstdlib>
#include<cstring>
using namespace std;

#define _DEBUG

#define MAXV 1005
#define MAXC 500
#define MAXE 85	//20*4+2//#define INF 1000000000 //此处建议不要过大或过小,过大易导致运算时溢出,过小可能会被判定为真正的距离
#define INF (1 << 15)
#define inputfile "case_example\\case1.txt"
#define outputfile "case_example\\ans1.txt"

const int CodeLen = MAXV;	 //定义编码的长度
const int GroupSize = 5;	//定义种群大小
const int Generations = 20;	//指定种群的繁殖代数 

struct edge;
struct point;
struct Chromosome;
//struct G;

extern vector<edge> G[];
extern vector<edge> Gout[];
extern vector<point> np;
extern int mincost;

//int spfa(vector<edge>(&G)[MAXV], int & V, int & sour, int *dist, int *prevv, int *preve, int *dprev, int *dpree, bool *done, int *consume);
int spfa(vector<edge>(&G)[MAXV], int (&dist)[MAXV], int (&prevv)[MAXV], int (&preve)[MAXV], bool (&done)[MAXV]);
//int min_cost_flow(vector<edge> (&G)[MAXV], int & sour, int & f, int & sink, int & currentcost, int *dist, int *prevv, int *preve, int *dprev, int *dpree, int *consume, bool *done, int & reflow);
//int min_cost_flow(vector<edge> (&G)[MAXV], int & f, int & currentcost, int (&dist)[MAXV], int (&prevv)[MAXV], int (&preve)[MAXV], bool (&done)[MAXV]);
int min_cost_flow(vector<edge> (&G)[MAXV], int & f, int & currentcost);
void InitialGroup(vector<Chromosome> & Group, int & V, int & C, int & needed, int & BC);
void Selection(vector<Chromosome> & Group);
void CrossOver(vector<Chromosome> & Group, int & V);
void Variation(vector<Chromosome> & Group, int & V);
bool comparebw(const point & np1, const point & np2);
void readdata(int & sour, int & sink, int & V, int & C, int & E, int & BC, int & needed, int (&consume)[MAXV], char * fin[MAX_EDGE_NUM], int line_num);
bool display(int & sour, int & sink, int *dprev, int *dpree, int * consume, int & totcost);

struct edge
{
	int to;	//连接的点
	int v;	//volume容量限制
	int c;	//cost带宽成本
	int rev;//reverse反向边
	edge(int to = 0, int volume = 0, int cost = 0, int reverse = 0) : to(to), v(volume), c(cost), rev(reverse) {};
};

struct point {
	int id;
	int bw, ct, du;//网络节点总带宽，总花费，度
	point(int id, int bandwidth = 0, int cost = 0, int du = 0) : id(id), bw(bandwidth), ct(cost), du(du) {};
};

struct Chromosome
{
	bool *key;    //用随机序列0/1表示安置/不安置服务器点
	int cost;
	int flow;	//未完成的流量，若为0即满足题目要求
	float fitvalue;      //适应值，用在选择遗传里面
	vector<edge> graph[MAXV];    //该染色体对应的邻接表

	Chromosome() { key = nullptr; flow = 0; cost = 0; fitvalue = 0.0; }
	Chromosome(int & V, int & ServeMax, int & needed, int & BC)
	{
		flow = needed;
		fitvalue = 0.0;
		for (int i = 0; i < V; i++) { graph[i].reserve(G[i].capacity()); graph[i] = G[i]; }

		do
		{
			key = new bool[V]();
			srand((unsigned)time(0));
			int ServeNum = rand() % ServeMax;  //服务器数量不超过消费节点数量
			for (; ServeNum > 0; ServeNum--) key[rand() % V] = 1;
			//zongdaikuan>zongxuqiu
			for (int t = 0; t < V; t++) {   
				if (key[t]) {
					cost += BC;   
					graph[V].push_back(edge(t, np[t].bw, 0, graph[t].size()));  //graph[V]是总源点
					graph[t].push_back(edge(V, 0, 0, graph[V].size() - 1));
				}
			}
			cost = INF;
			cost = min_cost_flow(graph, flow, cost);
			if (cost != -1) break;
			if (cost<mincost)
			{
				mincost = cost;
				for(int i=0;i<V;i++) Gout[i]=graph[i];
			}
			delete[]key;
		} while (true);
	}

	Chromosome(const Chromosome& o, int & V)
	{
		key = new bool[V]();
		for (int i = 0; i < V; i++) key[i] = o.key[i];
		cost = o.cost; flow = o.flow; fitvalue = o.fitvalue;
		for (int i = 0; i < V; i++) graph[i] = o.graph[i];
	}

	~Chromosome() { delete[]key; key = nullptr; }
	void print() { cout << key << endl; };	
};
///////////////////////////////////////////////////////
void deploy_server(char * graph[MAX_EDGE_NUM], int edge_num, char * filename);

	
#endif
