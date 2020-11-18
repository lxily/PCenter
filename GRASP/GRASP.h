#pragma once

#include <iostream>
#include <vector>
#include <ctime>
#include <map>
#include <cassert>
#include <algorithm>

#define inf 0x3f3f3f3f
#define pii pair<int,int>
#define piii pair<pii, int>
#define mkp(a, b) make_pair(a, b)

#define db1(a) cerr << #a << "=" << (a) << endl;
#define db2(a,b) cerr << #a << "=" << (a) << " " << #b << "=" << (b) << endl;
#define db3(a,b,c) cerr << #a << "=" << (a) <<" "<< #b << "=" << (b) << " " << #c << "=" << (c) << endl;

using namespace std;

using Solution = vector<int>;

enum {
	Greedy = false,
	Random = true,

	/*局部搜索参数*/
	TabuSearchSteps = 1000,
	NoDecreaseLimit = 1000,
	RandomAlphaInInit = 70,	// /100
	RandomAlphaInSearch = 85,	// /100

	/*PathRellinking参数*/
	DistBeta = 50,		// /100
	PopulationSize = 15,
	MaxIterations = 1000,

	/*禁忌表*/
	DoubleTabu = true,	//双禁忌效果更好
};

struct Rander {
	Rander(unsigned _seed = 0) {
		srand(_seed ? _seed : (unsigned)time(NULL));
	}
	int getRandNumber(int _l = -1, int _r = -1) {
		if (_l < 0 || _r < 0) {
			return rand();
		}
		if (_l < 0 || _r < _l) {
			cout << "RandomNumber Error!\n";
			exit(-1);
		}
		return rand() % (_r - _l + 1) + _l;
	}
};

struct Edge {
	int dist, from, to;
	Edge() {}
	Edge(int d, int f, int t) {
		dist = d; from = f; to = t;
	}
	void output() {
		cout << "(" << from << " -> " << to << ", dist=" << dist << ") ";
	}
	bool operator < (const Edge& E) const {
		return dist < E.dist || (dist == E.dist && to < E.to);
	}
};

struct TabuTable {
	int R, tableLen;
	bool doubleTabu;
	int vNum, pNum;		//For double Tabu
	vector<int>serviceTabu;
	vector<int>costomerTabu;
	vector<vector<int>>tabuTable;

	TabuTable(int _R = 10, int _vNum = 1, int _pNum = 1) {
		R = _R; tableLen = _vNum;
		doubleTabu = DoubleTabu;
		vNum = _vNum; pNum = _pNum;
		serviceTabu = vector<int>(tableLen, 0);		//成为服务节点的禁忌步长
		costomerTabu = vector<int>(tableLen, 0);		//...
		tabuTable = vector<vector<int>>(tableLen, vector<int>(tableLen, 0));
	}

	void clear() {
		serviceTabu = vector<int>(tableLen, 0);
		costomerTabu = vector<int>(tableLen, 0);
		tabuTable = vector<vector<int>>(tableLen, vector<int>(tableLen, 0));
	}

	bool isTabu(int u, int v, int nowStep) {	//u想加入成为服务节点
		if (doubleTabu) {
			return nowStep < serviceTabu[v] || nowStep < costomerTabu[u];
		}
		else {
			return nowStep < tabuTable[u][v];	//应该是不能变回<v, u>?
		}

	}

	void updateTabuTable(int u, int v, int nowStep, Rander rander) {	//u加入成为服务节点
		if (doubleTabu) {
			int alpha = nowStep + pNum / 10 + rander.getRandNumber() % 10;				//u alpha步之内不能再变回为普通节点
			int beta = nowStep + (vNum - pNum) / 10 + rander.getRandNumber() % 100;		//反之v Bate步之内不能成为服务节点
			costomerTabu[u] = alpha;
			serviceTabu[v] = beta;
		}
		else {
			//tabuTable[u][v] = nowStep + R + rander.getRandNumber(1, R+1);		//<u, v>被禁忌
			tabuTable[u][v] = nowStep + pNum * (vNum - pNum) / 100 + rander.getRandNumber(0, 10 * pNum);
		}
	}
};

struct PCenterSolver {

	///原始输入数据
	int vNum, pNum;	//vNum: 总节点数, pNum: 总中心数
	vector<vector<int>> graph;	//邻接矩阵表示的图

	///预处理数据
	vector<vector<Edge>> minDistEdge;	//minDistAndIndex[i][j] 表示距离点i出发第j长的边
	vector<vector<int>> kValue;			//kValue[i][j]=k 表示点j是距离点i第k近的点,距离相等取较小的小标

	///辅助数据结构
	vector<pii> FTable;
	vector<pii> DTable;
	TabuTable tabuTable;

	///计算结果
	int minMaximum;
	Solution centers;
	vector<bool>isCenter;

	int historyOptimalValue;
	Solution historyOptimalCenters;

	PCenterSolver() {};
	~PCenterSolver() {};
	PCenterSolver(vector<vector<int>> g, int p);	//初始化输入
	bool operator == (const PCenterSolver & sol)const;	//判断两个解是否相同(可优化)
	///算法过程

	void init(bool randomSolution = false);	//预处理
	Solution getInitSolution();	//得到一个初始解 O(p^2 * n)
	//返回距离当前解sol最远的点集合（所有的最长服务边）; 格式(u, v)->点u服务点v; 仅在得到初始解使用; 复杂度O(n*p)
	vector<pii> getMaxDistPointsFromSolution(const Solution &sol);
	vector<pii> getMaxDistPointsFromSolution();						//同上,区别在于通过查表得到，O(n)
	vector<pii> getNeighbourhoodsBySolution(const Solution &sol);	//返回当前解的交换邻域, O(N_k * p)	--已经舍弃
	vector<int> neighbourhoodsEvaluation(const vector<pii> &Neighbourhoods);	//邻域评估 O(N_k * p * n)	--已经舍弃
	vector<piii> getNeighbourhoodsAndEvaluation(const Solution &sol);	//邻域评估 O(N_k * n)
	int initFDTable(const Solution &sol);	//初始化FD表 O(n^2)
	int addPointToSolution(int addPoint);	//从当前解sol中增加一个点u，修改F表，返回Sc值 O(n)
	int removePointFromSolution(int removePoint); //从当前解sol中删除一个解u, 查询并修改D表，返回Mf值 O(n * p)
	int resultOfRemovePoint(int removePoint); //从当前解sol中删除一个解u, 查询D表，返回Mf值 O(n)
	void swapPairs(int addPoint, int removePoint);	//对当前解进行交换操作 O(n)
	void upDateHistoryOptimal();		//由自身解更新
	void upDateHistoryOptimal(const PCenterSolver &sol); //由指定解得到
	int searchOneStep(int nowStep);	//进行一步禁忌搜索 O(N_k * p * n + n^2) -> 邻域生成 + 评估
	void solveWithNSearch(int steps);	//进行N步禁忌搜索
	void changeToOptimalSolution();		//将历史最优解更新到当前解
	int calculateResultByForce(const Solution &sol);	//暴力计算结果，仅用于验算FDTable

	/// 输出 FOR DEBUG
	void printMinDistEdge();
	void printKValue();
	void printSolution(const Solution &sol);
	void printFTable();
	void printDTable();
};

struct PathRellinking {
	int a;
	struct QRecoder {
		int result, vc, vci, ve, vei;
	};
	PathRellinking() {}
	~PathRellinking() {}
	static vector<bool> pointsAppear(const vector<int>& points, int vNum); //O(n)
	static vector<int> intersectionOfSolution(const PCenterSolver& Sc, const PCenterSolver& Se); //O(n)
	static vector<int> removePointsFromSolution(const PCenterSolver& sol, const vector<int> &removePoints); //O(n)
	static vector<int> addPointsToSolution(const PCenterSolver& sol, const vector<int> &addPoints);  //O(n)
	static int distOfTwoSolution(const PCenterSolver& Sc, const PCenterSolver& Se);	//O(n)
	static bool solverInPopulation(const PCenterSolver& sol, const vector<PCenterSolver>& P);
	static PCenterSolver pathRellinkingCrossover(const PCenterSolver &Sc, const PCenterSolver &Se);
	static Solution GRASP(vector<vector<int>>_graph, int pNum, unsigned _randSeed = 0, int opt = 0);
};

extern Rander rander;
extern map<string, int>baselines;
extern vector<vector<string>>instPmed;