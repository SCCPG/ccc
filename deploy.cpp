#include "deploy.h"
#include <stdio.h>

vector<edge> G[MAXV];   //图的邻接表表示。由于一般情况下E<<V*V,故在此选用了vector动态数组存储,也可以使用链表存储
vector<edge> Gout[MAXV];
vector<point> np;	//网络节点
int mincost = INF;

using namespace std;
void InitialGroup(vector<Chromosome> & Group, int & V, int & C, int & needed, int & BC)
{
	/*


	用剪枝和搜索初始化种群


	*/
	
	/////////////////////////////////////////
	//不知道这里对内存的调用会不会出现问题///
	/////////////////////////////////////////
	for (int i = Group.size(); i < GroupSize; i++)
		Group.push_back(Chromosome(V, C, needed, BC));

	//打印初始种群
	for (vector<Chromosome>::iterator iter = Group.begin(); iter != Group.end(); ++iter) (*iter).print();
}

void Selection(vector<Chromosome> & Group)
{
	int CostSum = 0;
	float AccumuVal = 0.0, ptr = 0.0;
	vector<Chromosome> Orig(Group);
	Group.clear();

	for (vector<Chromosome>::iterator iter = Orig.begin(); iter != Orig.end(); iter++) 
		CostSum += (*iter).cost;
	for (vector<Chromosome>::iterator iter = Orig.begin(); iter != Orig.end(); iter++) 
		(*iter).fitvalue = (*iter).cost / (float)CostSum;

	srand((unsigned)time(0));
	for (int i = 0; i < GroupSize; i++) 
	{
		ptr = (rand() % 1000) / 1000.0;
		for (vector<Chromosome>::iterator iter = Orig.begin(); iter != Orig.end(); iter++) 
		{
			AccumuVal += (*iter).fitvalue;
			if (ptr <= AccumuVal) { Group.push_back((*iter)); break; }
		}
		AccumuVal = 0;
	}
	//打印选择遗传后的种群
	cout << "选择遗传" << endl;
	for (vector<Chromosome>::iterator iter = Group.begin(); iter != Group.end(); ++iter) 
		(*iter).print();
}

void CrossOver(vector<Chromosome> & Group, int & V)
{
	int rd, temp, sequence[GroupSize];
	srand((unsigned)time(0));

	for (int i = 0; i < GroupSize; i++) sequence[i] = i;
	for (int i = GroupSize - 1; i >= 0; i--)
	{
		rd = rand() % (i + 1);    // 从剩余的坐标中随机选取一个坐标，交换该坐标对应的值与i坐标对应的值 
		temp = sequence[rd];      // 这样得到的序列就是一个从0~N-1的随机打乱的序列
		sequence[rd] = sequence[i];
		sequence[i] = temp;
	}

	for (int i = 0; i < GroupSize - 1; i = i + 2)
	{
		int pos = rand() % V;    //选择要交叉的起始bit位
		cout << "序列" << sequence[i] << "和" << sequence[i + 1] << "的";
		cout << pos << "位交换." << endl;
		for (int j = pos; j < V; j++)
		{
			bool temp = Group.at(sequence[i]).key[j];
			Group.at(sequence[i]).key[j] = Group.at(sequence[i + 1]).key[j];
			(Group.at(sequence[i + 1]).key)[j] = temp;
		}
	}

	//打印交叉遗传后的种群
	cout << "交叉遗传" << endl;
	for (vector<Chromosome>::iterator iter = Group.begin(); iter != Group.end(); ++iter) (*iter).print();

}

void Variation(vector<Chromosome> & Group, int & V)
{
	int pos = 0;
	srand((unsigned)time(0));
	for (int i = 0; i < GroupSize; i++)
	{
		pos = rand() % V;  //选择要随机改变的bit位
		if (Group.at(i).key[pos]) Group.at(i).key[pos] = 0;
		else Group.at(i).key[pos] = 1;
	}

	//打印变异遗传后的种群
	cout << "变异遗传" << endl;
	for (vector<Chromosome>::iterator iter = Group.begin(); iter != Group.end(); ++iter) (*iter).print();
}

bool comparebw(const point & np1, const point & np2)
{
	return np1.bw > np2.bw;
}

//int spfa(vector<edge>(&G)[MAXV],int & V, int & sour, int *dist, int *prevv, int *preve, int *dprev, int *dpree, bool *done, int *consume)
int spfa(vector<edge>(&G)[MAXV], int (&dist)[MAXV], int (&prevv)[MAXV], int (&preve)[MAXV], bool (&done)[MAXV])
{ //返回值:TRUE为找到最短路返回,FALSE表示出现负环退出
    int V = G->size();
    int sink = V - 1;
    int sour = sink - 1;
	queue<int> buff;	
	/////////////////?????????????????/////////////////////
	//////           这里的要改成new  /////////////////////
	///////////////////////////////////////////////////////
	int cnt[MAXV];//记录入队次数,超过V则表示出现负环

	fill(dist, dist + V, INF);
	memset(cnt, 0, sizeof(cnt));
	memset(done, 0, sizeof(done));
	memset(prevv, -1, sizeof(prevv));
	dist[sour] = 0;
	buff.push(sour); //原点入队
	done[sour] = 1; //标记原点已经入队
	cnt[sour] = 1; //修改入队次数为1
	while (!buff.empty()){ //队列非空,需要继续松弛
		int tmp = buff.front(); //取出队首元素
		for (int i = 0; i<(int)G[tmp].size(); i++){ //枚举该边连接的每一条边
			edge *t = &G[tmp][i]; //由于vector的寻址速度较慢,故在此进行一次优化
			if ((*t).v>0 && (dist[tmp] + (*t).c<dist[(*t).to])){ //更改后距离更短,进行松弛操作
				dist[(*t).to] = dist[tmp] + (*t).c; //更改边权值
				prevv[(*t).to] = tmp;//路径还原//key
				preve[(*t).to] = i;//key
				if (!done[(*t).to]){ //没有入队,则将其入队
					buff.push((*t).to); //将节点压入队列
					done[(*t).to] = 1; //标记节点已经入队
					cnt[(*t).to] += 1; //节点入队次数自增
					if (cnt[(*t).to]>V){ //已经超过V次,出现负环
						while (!buff.empty())buff.pop(); //清空队列,释放内存
						return tmp; //返回FALSE
					}
				}
			}
		}
		buff.pop();//弹出队首节点
		done[tmp] = 0;//将队首节点标记为未入队
	}
	return -1; //返回-1
} //算法结束
//int spfa(vector<edge>(&G)[MAXV], int (&dist)[MAXV], int (&prevv)[MAXV], int (&preve)[MAXV], bool (&done)[MAXV])
//int min_cost_flow(vector<edge> (&G)[MAXV], int & f, int & currentcost, int (&dist)[MAXV], int (&prevv)[MAXV], int (&preve)[MAXV], bool (&done)[MAXV])
int min_cost_flow(vector<edge> (&G)[MAXV], int & f, int & currentcost)
{
    int V = G->size();
    int t = V - 1;
    int s = t - 1;
	int i;
	int dist[MAXV];
	int prevv[MAXV];
	int preve[MAXV];
	bool done[MAXV];
	//currentcost = 0;
	while (f > 0)
	{
		int d = f; //d:本次求得的最小费用流
		int p = spfa(G,dist,prevv,preve,done);
		if (p == -1)
		{
			if (dist[t] == INF){
				//reflow = f;
				return -1;
			}
			for (i = t; i != s; i = prevv[i])
				d = min(d, G[prevv[i]][preve[i]].v);

			currentcost += d * dist[t];
			f -= d;
			//cout << f << '\n';
			for (i = t; i != s; i = prevv[i])
			{
				edge &es = G[prevv[i]][preve[i]];
				es.v -= d;
				G[es.to][es.rev].v += d;
			}
		}
		else
		{
			memset(done, 0, sizeof(done));
			while (!done[p]) {
				done[p] = 1;
				p = prevv[p];
			}
			int c = G[prevv[p]][preve[p]].c;
			d = min(d, G[prevv[p]][preve[p]].v);
			for (i = prevv[p]; i != p; i = prevv[i]){
				d = min(d, G[prevv[i]][preve[i]].v);
				c += G[prevv[i]][preve[i]].c;
			}
			if (d == 0)
				return currentcost;
			currentcost += d * c;
			edge &es = G[prevv[p]][preve[p]];
			es.v -= d;
			G[es.to][es.rev].v += d;
			for (i = prevv[p]; i != p; i = prevv[i])
			{
				edge &es = G[prevv[i]][preve[i]];
				es.v -= d;
				G[es.to][es.rev].v += d;
			}
		}
	}
	return currentcost;
}

void readdata(int & sour, int & sink, int & V, int & C, int & E, int & BC, int & needed, int (&consume)[MAXV], char * fin[MAX_EDGE_NUM], int line_num)//读入网络节点和消费节点，消费节点用总汇点sink代替，总源点未连接
{
	int cnt=0;
	sscanf(fin[cnt++], "%d%d%d", &V, &E, &C); //读入点数和边数
	sour = V++;
	sink = V++;
	sscanf(fin[cnt++], "%d", &BC);	//读入部署成本
	np.reserve(sour);
	for (int i = 0; i < sour; i++){ G[i].reserve(MAXE); np.push_back(point(i)); }
	G[sour].reserve(V); G[sink].reserve(V);
	for (int i = 0, x, y, l, c; i < E; i++){
		sscanf(fin[cnt++], "%d%d%d%d", &x, &y, &l, &c); //读入x,y,l表示从x->y有一条有带宽限制为l，带宽成本为c的边
		G[x].push_back(edge(y, l, c, G[y].size()));
		G[y].push_back(edge(x, 0, -c, G[x].size() - 1));
		G[x].push_back(edge(y, 0, -c, G[y].size()));
		G[y].push_back(edge(x, l, c, G[x].size() - 1));
		np[x].bw += l;	np[x].ct += c*l;	np[x].du++;
		np[y].bw += l;	np[y].ct += c*l;	np[y].du++;
	}
	needed = 0;
	for (int i = 0, x, y, vol; i < C; i++){
		sscanf(fin[cnt++], "%d%d%d", &x, &y, &vol); //读入x,y,l表示从x->y有一条有带宽限制为l，带宽成本为c的边
		consume[y] = x;
		needed += vol;
		G[y].push_back(edge(sink, vol, 0, G[sink].size()));
		G[sink].push_back(edge(y, 0, 0, G[y].size() - 1));
		np[y].bw += vol;	np[y].du++;
	}
}
/*
bool display(int & sour, int & sink, int *dprev, int *dpree, int * consume, int & totcost)
{
	FILE *fout;
	if ((fout = fopen(outputfile, "w")) != 0){
		printf("The file spfax.txt was not opened\n");
		return false;
	}
	fprintf(fout, "                  \n");
	int paths = 0, minvol = INF, cost = 0, c = 0;
	unsigned tmp = sour, tmpe = 0, tmpv = 0;
	cost = 0;
	for (tmpe = 0;;)
	{
		if ((tmp == sour) && (tmpe >= Gout[tmp].size()))
			break;
		edge& e = Gout[tmp][tmpe];
		if (e.c<0){
			tmpe++;
			continue;
		}
		if (Gout[e.to][e.rev].v > 0)
		{
			cost += e.c;
			tmpv = e.to;
			minvol = min(minvol, Gout[e.to][e.rev].v);
			if (tmpv != sink){
				fprintf(fout, "%d\t", tmpv);
				dprev[tmpv] = tmp;
				dpree[tmpv] = tmpe;
				tmp = tmpv;
				tmpe = 0;
			}
			else
			{
				c = cost*minvol;
				totcost += c;
				fprintf(fout, "%d\t%d\t%d\n", consume[tmp], minvol, c);
				//Gout[e.to][e.rev].v -= minvol;
				dprev[tmpv] = tmp;
				dpree[tmpv] = tmpe;
				for (int i = tmpv; i != sour; i = dprev[i])
				{
					edge &es = Gout[dprev[i]][dpree[i]];
					//es.v -= minvol;
					Gout[es.to][es.rev].v -= minvol;
				}
				tmp = sour;
				tmpe = 0; 
				minvol = INF;
				cost = 0;
				paths++;
			}
		}
		else
			tmpe++;
	}
	fprintf(fout, "\n%d", totcost);
	rewind(fout);
	fprintf(fout, "%d\n", paths);
	fclose(fout);
	return true;
}
*/
//你要完成的功能总入口
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
	int V;	//网络节点数
	int E;	//边数
	int C;	//消费节点数
	int BC;	//部署成本
	int sour;
	int sink;
	int needed;
	int reflow;
	int totcost;
	int currentcost;

	int dist[MAXV]; //存储到原点0的距离,可以开二维数组存储每对节点之间的距离
	int prevv[MAXV], preve[MAXV];//最短路中的前驱结点 对应边
	int consume[MAXV];//网络节点连接的消费节点编号
	int dprev[MAXV], dpree[MAXV];//输出路径中的前驱结点 对应边
	bool done[MAXV]; //用于判断该节点是否已经在队列中
	queue<int> buff; //队列,用于存储在SPFA算法中的需要松弛的节点

	vector<Chromosome> Group;	//遗传种群

	Group.reserve(GroupSize);
	readdata(sour,sink,V,C,E,BC,needed,consume,topo,line_num);
	InitialGroup(Group,V,C,needed,BC);
	//Selection(Group);
	//CrossOver(Group,V);
	//Variation(Group,V);

	// 需要输出的内容
	char * topo_file = (char *)"17\n\n0 8 0 20\n21 8 0 20\n9 11 1 13\n21 22 2 20\n23 22 2 8\n1 3 3 11\n24 3 3 17\n27 3 3 26\n24 3 3 10\n18 17 4 11\n1 19 5 26\n1 16 6 15\n15 13 7 13\n4 5 8 18\n2 25 9 15\n0 7 10 10\n23 24 11 23";

	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
	write_result(topo_file, filename);

}
