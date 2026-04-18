// Taken from https://github.com/ekzhang/library/blob/master/maxflow2.cpp and modified
// maxflow2.cpp
// Eric K. Zhang; Aug. 7, 2017

using namespace std;

/* Maximum-Flow solver using Dinic's Blocking Flow Algorithm.
 * Time Complexity:
 *   - O(V^2 E) for general graphs, but in practice ~O(E^1.5)
 *   - O(V^(1/2) E) for bipartite matching
 *   - O(min(V^(2/3), E^(1/2)) E) for unit capacity graphs
 */
template<class T>
class max_flow {
	static const T INF = numeric_limits<T>::max();

	struct edge {
		int t, rev;
		T cap, f;
	};

	std::vector<std::vector<edge>> adj;
	std::vector<int> dist;
	std::vector<int> ptr;
	int m_size;

	bool bfs(int s, int t) {
		std::fill(dist.begin(), dist.end(), -1);
		dist[s] = 0;
		queue<int> q({ s });
		while (!q.empty() && dist[t] == -1) {
			int n = q.front();
			q.pop();
			for (auto& e : adj[n]) {
				if (dist[e.t] == -1 && e.cap != e.f) {
					dist[e.t] = dist[n] + 1;
					q.push(e.t);
				}
			}
		}
		return dist[t] != -1;
	}

	T augment(int n, T amt, int t) {
		if (n == t) return amt;
		for (; ptr[n] < (int)adj[n].size(); ptr[n]++) {
			edge& e = adj[n][ptr[n]];
			if (dist[e.t] == dist[n] + 1 && e.cap != e.f) {
				T flow = augment(e.t, min(amt, e.cap - e.f), t);
				if (flow != 0) {
					e.f += flow;
					adj[e.t][e.rev].f -= flow;
					return flow;
				}
			}
		}
		return 0;
	}

public:
	max_flow(int size): m_size(size){
		adj.resize(m_size);
		dist.resize(m_size);
		ptr.resize(m_size);
	}

	void add(int u, int v, T cap=1, T rcap=0) {
		if (u >= m_size || v >= m_size) {
			throw std::out_of_range("Node index out of range in add()");
		}
		adj[u].push_back({ v, (int) adj[v].size(), cap, 0 });
		adj[v].push_back({ u, (int) adj[u].size() - 1, rcap, 0 });
	}

	T calc(int s, int t) {
		T flow = 0;
		while (bfs(s, t)) {
			std::fill(ptr.begin(), ptr.end(), 0);
			while (T df = augment(s, INF, t))
				flow += df;
		}
		return flow;
	}

	void clear() {
		for (int n = 0; n < m_size; n++)
			adj[n].clear();
	}
};
