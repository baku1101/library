template<typename T>
class bellmanford {
public:
	struct edge {
		int from, to;
		T cost;
	};
	int V;
	T inf;
	vector<T> d;
	vector<edge> es;
	bellmanford(int n) : V(n), inf(numeric_limits<T>::max()/2), d(n, inf){}
	void add_edge(int from, int to, T cost){
		es.push_back((edge){from,to,cost});
	}
	bool solve(int s) {
		d[s] = 0;
		T cnt;
		for (cnt = 0; cnt < V; cnt++) {
			bool update = false;
			for (auto&& e : es) {
				if (d[e.from] < inf and d[e.to]  > d[e.from] + e.cost) {
					d[e.to] = d[e.from] + e.cost;
					update = true;
				}
			}
			if (!update) {
				break;
			}
		}
		return (cnt == V);
	}
};
