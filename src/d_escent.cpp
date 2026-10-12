/*
 * D-ESCENT: Dynamic Edge-Independent Spanning Trees
 *
 * Compile: g++ -O3 -fopenmp -std=c++17 -o d_escent src/d_escent.cpp
 * Usage:   ./d_escent <dataset_file> [--validate]
 */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// ===== Timer =====
struct Timer {
    using Clock = chrono::high_resolution_clock;
    Clock::time_point s;
    void start() { s = Clock::now(); }
    double ms() const { return chrono::duration<double,milli>(Clock::now()-s).count(); }
};

// ===== Canonical edge =====
struct Edge {
    int u, v;
    Edge() : u(-1), v(-1) {}
    Edge(int a, int b) : u(min(a,b)), v(max(a,b)) {}
    bool operator==(const Edge& o) const { return u==o.u && v==o.v; }
    bool valid() const { return u >= 0; }
};
struct EdgeHash {
    size_t operator()(const Edge& e) const {
        return hash<long long>()((long long)e.u * 1000003LL + e.v);
    }
};
using EdgeSet = unordered_set<Edge, EdgeHash>;

// ===== Simple Graph (adjacency lists + edge set) =====
struct Graph {
    int N = 0;
    vector<vector<int>> adj;
    vector<bool> active;
    EdgeSet edges;

    void init(int n) {
        N = n; adj.assign(n, {}); active.assign(n, false);
    }
    void grow(int n) {
        if (n <= N) return;
        adj.resize(n); active.resize(n, false); N = n;
    }
    void add_edge(int u, int v) {
        grow(max(u,v)+1);
        active[u] = active[v] = true;
        Edge e(u,v);
        if (edges.insert(e).second) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
    void remove_edge(int u, int v) {
        Edge e(u,v);
        if (edges.erase(e)) {
            adj[u].erase(remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
            adj[v].erase(remove(adj[v].begin(), adj[v].end(), u), adj[v].end());
        }
    }
    bool has_edge(int u, int v) const { return edges.count(Edge(u,v)); }
    int num_verts() const { int c=0; for(int i=0;i<N;i++) if(active[i]) c++; return c; }
    int num_edges() const { return (int)edges.size(); }
    vector<int> verts() const { vector<int> r; for(int i=0;i<N;i++) if(active[i]) r.push_back(i); return r; }
};

// ===== Tarjan BiCC: extract largest biconnected component =====
// Returns vertex set of the largest BiCC in the graph.
vector<int> largest_bicc_verts(const vector<vector<int>>& adj, const vector<bool>& active, int N) {
    vector<int> disc(N,-1), low(N,-1), par(N,-1);
    int timer_v = 0;

    // Store components as list of vertex-sets
    vector<vector<int>> comps;

    // Edge stack for BiCC extraction
    struct SE { int u, v; };
    stack<SE> estk;

    auto pop_comp = [&](int u, int v) {
        unordered_set<int> vs;
        while (!estk.empty()) {
            auto e = estk.top(); estk.pop();
            vs.insert(e.u); vs.insert(e.v);
            if ((e.u==u && e.v==v) || (e.u==v && e.v==u)) break;
        }
        if (vs.size() >= 2)
            comps.push_back(vector<int>(vs.begin(), vs.end()));
    };

    // Iterative DFS with explicit stack
    struct Frame { int u; int idx; bool returning; int child; };

    for (int start = 0; start < N; start++) {
        if (!active[start] || disc[start] != -1) continue;

        stack<Frame> stk;
        disc[start] = low[start] = timer_v++;
        stk.push({start, 0, false, -1});

        while (!stk.empty()) {
            auto& fr = stk.top();
            int u = fr.u;

            if (fr.returning) {
                int v = fr.child;
                low[u] = min(low[u], low[v]);
                if (low[v] >= disc[u]) {
                    pop_comp(u, v);
                }
                fr.returning = false;
                fr.idx++;
            }

            bool pushed = false;
            while (fr.idx < (int)adj[u].size()) {
                int v = adj[u][fr.idx];
                if (disc[v] == -1) {
                    par[v] = u;
                    disc[v] = low[v] = timer_v++;
                    estk.push({u, v});
                    fr.child = v;
                    fr.returning = true;
                    stk.push({v, 0, false, -1});
                    pushed = true;
                    break;
                } else if (v != par[u] && disc[v] < disc[u]) {
                    estk.push({u, v});
                    low[u] = min(low[u], disc[v]);
                }
                fr.idx++;
            }

            if (!pushed && !fr.returning) {
                stk.pop();
            }
        }

        // Remaining edges form a component
        if (!estk.empty()) {
            unordered_set<int> vs;
            while (!estk.empty()) {
                auto e = estk.top(); estk.pop();
                vs.insert(e.u); vs.insert(e.v);
            }
            if (vs.size() >= 2)
                comps.push_back(vector<int>(vs.begin(), vs.end()));
        }
    }

    // Return largest by vertex count
    int best = -1, bestSz = 0;
    for (int i = 0; i < (int)comps.size(); i++) {
        if ((int)comps[i].size() > bestSz) {
            bestSz = (int)comps[i].size();
            best = i;
        }
    }
    if (best < 0) return {};
    return comps[best];
}

Graph extract_largest_bicc(const Graph& G) {
    auto verts = largest_bicc_verts(G.adj, G.active, G.N);
    Graph H;
    H.init(G.N);
    unordered_set<int> inV(verts.begin(), verts.end());
    for (int v : verts) H.active[v] = true;
    for (auto& e : G.edges) {
        if (inV.count(e.u) && inV.count(e.v))
            H.add_edge(e.u, e.v);
    }
    return H;
}

// ===== Even-Tarjan st-numbering (OGDF pathfinder approach) =====
// Precondition: H biconnected, s and t adjacent.
// Returns g[] with g[s]=1, g[t]=n, g[v] in {1..n} for all v in H.
vector<int> compute_st_numbering(const Graph& H, int s, int t) {
    int N = H.N;
    auto verts = H.verts();
    int n = (int)verts.size();
    vector<int> g(N, 0);

    if (n <= 0) return g;
    if (n == 1) { g[s] = 1; return g; }
    if (n == 2) { g[s] = 1; g[t] = 2; return g; }

    // Step 1: DFS. Pre-number t, then DFS from s.
    vector<int> dfn(N, 0);
    vector<int> low(N, 0);
    vector<int> dfsParent(N, -1);
    vector<int> followLow(N, -1);
    int cnt = 1;

    dfn[t] = cnt++; low[t] = dfn[t];
    dfn[s] = cnt++; low[s] = dfn[s];

    {
        struct F { int u; int idx; };
        stack<F> stk;
        stk.push({s, 0});
        while (!stk.empty()) {
            auto& fr = stk.top();
            int u = fr.u;
            bool pushed = false;
            while (fr.idx < (int)H.adj[u].size()) {
                int w = H.adj[u][fr.idx]; fr.idx++;
                if (dfn[w] == 0) {
                    dfsParent[w] = u;
                    dfn[w] = cnt++; low[w] = dfn[w];
                    stk.push({w, 0});
                    pushed = true;
                    break;
                } else if (low[u] > dfn[w]) {
                    low[u] = dfn[w];
                    followLow[u] = w;
                }
            }
            if (!pushed) {
                stk.pop();
                int p = dfsParent[u];
                if (p >= 0 && low[p] > low[u]) {
                    low[p] = low[u];
                    followLow[p] = u;
                }
            }
        }
    }
    if (low[t] > low[s]) low[t] = low[s];

    // Step 2: Pathfinder
    vector<bool> markedNode(N, false);
    markedNode[s] = markedNode[t] = true;

    // Edge marking
    EdgeSet markedEdgeSet;
    markedEdgeSet.insert(Edge(s, t));

    vector<int> adjIdx(N, 0);
    stack<int> nodeStack;
    nodeStack.push(t);
    nodeStack.push(s);

    int num = 1;
    int v = nodeStack.top(); nodeStack.pop();

    while (v != t) {
        bool pathFound = false;
        while (adjIdx[v] < (int)H.adj[v].size()) {
            int w = H.adj[v][adjIdx[v]]; adjIdx[v]++;
            Edge ew(v, w);
            if (markedEdgeSet.count(ew)) continue;
            markedEdgeSet.insert(ew);

            if (dfsParent[w] == v) {
                vector<int> path;
                path.push_back(v);
                while (!markedNode[w]) {
                    path.push_back(w);
                    markedNode[w] = true;
                    int nxt = followLow[w];
                    if (nxt < 0) break;
                    markedEdgeSet.insert(Edge(w, nxt));
                    w = nxt;
                }
                for (int i = (int)path.size()-1; i >= 0; i--)
                    nodeStack.push(path[i]);
                pathFound = true;
                break;
            } else if (dfn[v] < dfn[w]) {
                vector<int> path;
                path.push_back(v);
                while (!markedNode[w]) {
                    path.push_back(w);
                    markedNode[w] = true;
                    int par = dfsParent[w];
                    if (par < 0) break;
                    markedEdgeSet.insert(Edge(w, par));
                    w = par;
                }
                for (int i = (int)path.size()-1; i >= 0; i--)
                    nodeStack.push(path[i]);
                pathFound = true;
                break;
            }
        }
        if (!pathFound) {
            g[v] = num++;
            adjIdx[v] = 0;
        }
        v = nodeStack.top(); nodeStack.pop();
    }
    g[t] = num++;

    return g;
}

// ===== Itai-Rodeh tree construction =====
struct TreePair {
    vector<int> par1, par2; // parent arrays, -1 for root
    EdgeSet edg1, edg2;
    int root, t_vert;
};

TreePair build_trees(const Graph& H, const vector<int>& g, int root, int t_vert) {
    int N = H.N;
    TreePair tp;
    tp.root = root;
    tp.t_vert = t_vert;
    tp.par1.assign(N, -1);
    tp.par2.assign(N, -1);

    auto verts = H.verts();

    // T1: low-parent tree (parent has largest g < g[v])
    for (int v : verts) {
        if (v == root) continue;
        int best = -1, bestG = -1;
        for (int u : H.adj[v]) {
            if (g[u] > 0 && g[u] < g[v] && g[u] > bestG) {
                bestG = g[u]; best = u;
            }
        }
        tp.par1[v] = best;
        if (best >= 0) tp.edg1.insert(Edge(v, best));
    }

    // T2: high-parent tree (parent has smallest g > g[v])
    // Special: t's parent is root
    for (int v : verts) {
        if (v == root) continue;
        if (v == t_vert) {
            tp.par2[v] = root;
            tp.edg2.insert(Edge(v, root));
            continue;
        }
        int best = -1, bestG = INT_MAX;
        for (int u : H.adj[v]) {
            if (g[u] > 0 && g[u] > g[v] && g[u] < bestG) {
                bestG = g[u]; best = u;
            }
        }
        tp.par2[v] = best;
        if (best >= 0) tp.edg2.insert(Edge(v, best));
    }

    return tp;
}

// ===== BFS utilities =====
int bfs_farthest(const Graph& H, int start) {
    vector<int> dist(H.N, -1);
    queue<int> q;
    dist[start] = 0; q.push(start);
    int far = start, mx = 0;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : H.adj[u]) {
            if (dist[v]<0) {
                dist[v] = dist[u]+1;
                q.push(v);
                if (dist[v]>mx) { mx=dist[v]; far=v; }
            }
        }
    }
    return far;
}

struct TreeStats { int height; double avg; };
TreeStats tree_stats(const vector<int>& parent, int root, const vector<int>& verts, int N) {
    vector<vector<int>> ch(N);
    for (int v : verts) if (v!=root && parent[v]>=0) ch[parent[v]].push_back(v);
    vector<int> depth(N, -1);
    queue<int> q;
    depth[root]=0; q.push(root);
    int mx=0; long long tot=0; int cnt=0;
    while (!q.empty()) {
        int u=q.front(); q.pop(); cnt++; tot+=depth[u];
        if (depth[u]>mx) mx=depth[u];
        for (int c : ch[u]) { depth[c]=depth[u]+1; q.push(c); }
    }
    return {mx, cnt>0?(double)tot/cnt:0.0};
}

double compute_cost(const TreePair& tp, const vector<int>& verts, int N) {
    auto s1 = tree_stats(tp.par1, tp.root, verts, N);
    auto s2 = tree_stats(tp.par2, tp.root, verts, N);
    return max(s1.height,s2.height) + s1.avg + s2.avg;
}

int double_sweep_bfs(const Graph& H, int K) {
    auto verts = H.verts();
    if (verts.empty()) return -1;
    if (verts.size()==1) return verts[0];

    mt19937 rng(42);
    int trials = min(K,(int)verts.size());
    vector<int> starts(trials);
    for (int iter = 0; iter < trials; iter++)
        starts[iter] = verts[rng()%verts.size()];

    int bestRoot = verts[0];
    double bestCost = 1e18;

    #pragma omp parallel
    {
    int localRoot = bestRoot;
    double localCost = 1e18;

    #pragma omp for schedule(dynamic)
    for (int iter = 0; iter < trials; iter++) {
        int start = starts[iter];
        int v1 = bfs_farthest(H, start);
        // v1 as candidate
        if (H.adj[v1].empty()) continue;
        int tc = H.adj[v1][0];
        auto gc = compute_st_numbering(H, v1, tc);

        // Check validity
        bool ok = true;
        for (int v : verts) if (gc[v]<1||gc[v]>(int)verts.size()) { ok=false; break; }
        if (!ok) continue;

        auto tpc = build_trees(H, gc, v1, tc);
        bool valid = true;
        for (int v : verts) {
            if (v==v1) continue;
            if (tpc.par1[v]<0 || tpc.par2[v]<0) { valid=false; break; }
        }
        if (!valid) continue;

        double cost = compute_cost(tpc, verts, H.N);
        if (cost < localCost) { localCost = cost; localRoot = v1; }
    }

    #pragma omp critical
    {
        if (localCost < bestCost) {
            bestCost = localCost;
            bestRoot = localRoot;
        }
    }
    }
    return bestRoot;
}

// ===== D-ESCENT system =====
struct D_ESCENT {
    Graph G, H;
    TreePair tp;
    vector<int> g;
    int root=-1, t_vert=-1;

    // BiCC label per vertex
    vector<int> bicc_lbl;
    int next_lbl = 0;
    int main_lbl = -1;

    // Params
    int TAU = 50;
    int MERGE_TAU = 512;
    double alpha = 1.5;
    double growth_thr = 2.0;

    // Stats
    int n_ins=0, n_del=0, n_rebuild=0, n_fast=0, n_slow=0;

    EdgeSet non_tree;
    vector<EdgeSet> incident_non_tree;

    struct Replacement {
        Edge edge;
        int child = -1;
        int parent = -1;
        bool valid() const { return edge.valid() && child >= 0 && parent >= 0; }
    };

    // DS4
    double cost_last = 0;
    int upd_ctr=0, merge_ctr=0, size_last=0;
    int g_size_last=0;
    bool pending_rebuild=false;

    void init(int maxV) {
        G.init(maxV);
        bicc_lbl.assign(maxV, -1);
        incident_non_tree.assign(maxV, {});
    }

    void ensure_vertex_capacity(int n) {
        if ((int)bicc_lbl.size() < n) bicc_lbl.resize(n, -1);
        if ((int)incident_non_tree.size() < n) incident_non_tree.resize(n);
    }

    void add_non_tree_edge(const Edge& e) {
        non_tree.insert(e);
        ensure_vertex_capacity(max(e.u, e.v) + 1);
        incident_non_tree[e.u].insert(e);
        incident_non_tree[e.v].insert(e);
    }

    void remove_non_tree_edge(const Edge& e) {
        non_tree.erase(e);
        if (e.u < (int)incident_non_tree.size()) incident_non_tree[e.u].erase(e);
        if (e.v < (int)incident_non_tree.size()) incident_non_tree[e.v].erase(e);
    }

    void build_initial(const vector<pair<int,int>>& el) {
        for (auto& [u,v] : el) G.add_edge(u,v);
        do_full_rebuild(32);
    }

    void do_full_rebuild(int K) {
        H = extract_largest_bicc(G);
        auto verts = H.verts();
        if ((int)verts.size()<2) { cerr<<"BiCC < 2 verts\n"; return; }
        n_rebuild++;

        main_lbl = next_lbl++;
        for (int v : verts) bicc_lbl[v] = main_lbl;

        root = double_sweep_bfs(H, K);
        if (H.adj[root].empty()) root = verts[0];
        t_vert = H.adj[root][0];

        g = compute_st_numbering(H, root, t_vert);
        tp = build_trees(H, g, root, t_vert);

        non_tree.clear();
        incident_non_tree.assign(G.N, {});
        for (auto& e : H.edges) {
            if (!tp.edg1.count(e) && !tp.edg2.count(e))
                add_non_tree_edge(e);
        }

        cost_last = compute_cost(tp, verts, H.N);
        upd_ctr = 0; merge_ctr = 0;
        size_last = (int)verts.size();
        g_size_last = G.num_verts();
        pending_rebuild = false;
    }

    void rebuild_H(int K, int pref_root=-1) {
        Graph newH = extract_largest_bicc(H);
        auto verts = newH.verts();
        if ((int)verts.size()<2) return;
        H = newH;
        n_rebuild++;

        if (pref_root>=0 && pref_root<H.N && H.active[pref_root])
            root = pref_root;
        else
            root = double_sweep_bfs(H, K);
        if (H.adj[root].empty()) root = verts[0];
        t_vert = H.adj[root][0];

        g = compute_st_numbering(H, root, t_vert);
        tp = build_trees(H, g, root, t_vert);

        non_tree.clear();
        incident_non_tree.assign(G.N, {});
        for (auto& e : H.edges)
            if (!tp.edg1.count(e) && !tp.edg2.count(e))
                add_non_tree_edge(e);

        main_lbl = next_lbl++;
        for (int v : verts) bicc_lbl[v] = main_lbl;

        cost_last = compute_cost(tp, verts, H.N);
        size_last = (int)verts.size();
        g_size_last = G.num_verts();
        pending_rebuild = false;
        merge_ctr = 0;
    }

    bool same_bicc(int u, int v) {
        if (u<0||v<0||u>=(int)bicc_lbl.size()||v>=(int)bicc_lbl.size()) return false;
        return bicc_lbl[u]>=0 && bicc_lbl[u]==bicc_lbl[v];
    }

    Replacement find_repl(int child, bool isT1, Edge excl=Edge()) {
        if (child < 0 || child >= (int)incident_non_tree.size()) return {};
        for (const auto& e : incident_non_tree[child]) {
            if (excl.valid() && e == excl) continue;
            int other = (e.u == child) ? e.v : e.u;
            if (other < 0 || other >= H.N) continue;
            if (g[other] == 0 || g[child] == 0) continue;

            if (isT1 && g[other] < g[child]) return {e, child, other};
            if (!isT1 && g[other] > g[child]) return {e, child, other};
        }
        return {};
    }

    void insert_edge(int u, int v) {
        ensure_vertex_capacity(max(u,v)+1);
        Edge e(u,v);
        if (G.has_edge(u,v)) return;
        n_ins++; G.add_edge(u,v); upd_ctr++;

        if (same_bicc(u,v)) {
            H.add_edge(u,v);
            if (!tp.edg1.count(e) && !tp.edg2.count(e))
                add_non_tree_edge(e);
            n_fast++;
        } else {
            merge_ctr++;
            bool lazy = (merge_ctr<MERGE_TAU) && (G.num_verts()<(int)(growth_thr*g_size_last));
            if (lazy) { pending_rebuild=true; n_fast++; }
            else { do_full_rebuild(5); n_slow++; }
        }
        periodic();
    }

    void delete_edge(int u, int v) {
        n_del++;
        if (!G.has_edge(u,v)) return;
        G.remove_edge(u,v); upd_ctr++;
        Edge e(u,v);

        if (!H.has_edge(u,v)) { n_fast++; periodic(); return; }
        H.remove_edge(u,v);

        bool i1 = tp.edg1.count(e), i2 = tp.edg2.count(e);

        if (!i1 && !i2) {
            remove_non_tree_edge(e);
            n_fast++;
        } else if (i1 && !i2) {
            tp.edg1.erase(e);
            int child = -1;
            if (tp.par1[u]==v) child=u;
            else if (tp.par1[v]==u) child=v;
            else child = (g[u]>g[v])?u:v;
            Replacement r = find_repl(child, true);
            if (r.valid()) {
                tp.par1[r.child] = r.parent;
                tp.edg1.insert(r.edge);
                remove_non_tree_edge(r.edge);
                n_fast++;
            } else { rebuild_H(1,root); n_slow++; }
        } else if (!i1 && i2) {
            tp.edg2.erase(e);
            int child = -1;
            if (tp.par2[u]==v) child=u;
            else if (tp.par2[v]==u) child=v;
            else child = (g[u]<g[v])?u:v;
            Replacement r = find_repl(child, false);
            if (r.valid()) {
                tp.par2[r.child] = r.parent;
                tp.edg2.insert(r.edge);
                remove_non_tree_edge(r.edge);
                n_fast++;
            } else { rebuild_H(1,root); n_slow++; }
        } else {
            tp.edg1.erase(e); tp.edg2.erase(e);
            int c1=-1, c2=-1;
            if (tp.par1[u]==v) c1=u; else if (tp.par1[v]==u) c1=v;
            else c1 = (g[u]>g[v])?u:v;
            if (tp.par2[u]==v) c2=u; else if (tp.par2[v]==u) c2=v;
            else c2 = (g[u]<g[v])?u:v;
            Replacement r1 = find_repl(c1, true);
            Replacement r2 = find_repl(c2, false, r1.edge);
            if (r1.valid() && r2.valid()) {
                tp.par1[r1.child] = r1.parent;
                tp.edg1.insert(r1.edge); remove_non_tree_edge(r1.edge);
                tp.par2[r2.child] = r2.parent;
                tp.edg2.insert(r2.edge); remove_non_tree_edge(r2.edge);
                n_fast++;
            } else { rebuild_H(1,root); n_slow++; }
        }
        periodic();
    }

    void periodic() {
        if (upd_ctr < TAU) return;
        upd_ctr = 0;
        if (pending_rebuild) {
            bool merge_budget_full = merge_ctr >= MERGE_TAU;
            bool graph_grew = G.num_verts() >= (int)(growth_thr * g_size_last);
            if (merge_budget_full || graph_grew) {
                do_full_rebuild(5);
            }
            return;
        }
        auto verts = H.verts();
        if ((int)verts.size()>=2) {
            double cc = compute_cost(tp, verts, H.N);
            if (cc > alpha*cost_last) rebuild_H(5);
            else cost_last = cc;
        }
    }

    bool validate() {
        auto verts = H.verts();
        int n = (int)verts.size();
        if (n<2) return true;
        for (int v : verts) {
            if (v==root) continue;
            if (tp.par1[v]<0) { cerr<<"FAIL: T1 no parent v="<<v<<"\n"; return false; }
            if (g[tp.par1[v]]>=g[v]) { cerr<<"FAIL: T1 low-parent v="<<v<<"\n"; return false; }
            if (tp.par2[v]<0) { cerr<<"FAIL: T2 no parent v="<<v<<"\n"; return false; }
            if (v != t_vert && g[tp.par2[v]]<=g[v]) { cerr<<"FAIL: T2 high-parent v="<<v<<" g="<<g[v]<<" par="<<tp.par2[v]<<" g_par="<<g[tp.par2[v]]<<"\n"; return false; }
            if (v == t_vert && tp.par2[v]!=root) { cerr<<"FAIL: T2 t_vert parent not root\n"; return false; }
        }
        // Independence check on sample
        mt19937 rng(42);
        auto sv = verts;
        shuffle(sv.begin(), sv.end(), rng);
        int chk = min(n, 200);
        for (int i = 0; i < chk; i++) {
            int v = sv[i];
            if (v==root) continue;
            EdgeSet p1, p2;
            for (int c=v,d=0; c!=root&&d<n; d++) {
                int p=tp.par1[c]; if(p<0) break;
                p1.insert(Edge(c,p)); c=p;
            }
            for (int c=v,d=0; c!=root&&d<n; d++) {
                int p=tp.par2[c]; if(p<0) break;
                p2.insert(Edge(c,p)); c=p;
            }
            for (auto& e : p1) if (p2.count(e)) {
                cerr<<"FAIL: Independence v="<<v<<" edge("<<e.u<<","<<e.v<<")\n";
                return false;
            }
        }
        return true;
    }

    void print_stats() {
        cout<<"  |V(H)|="<<H.num_verts()<<" |E(H)|="<<H.num_edges()
            <<" |V(G)|="<<G.num_verts()<<" |E(G)|="<<G.num_edges()<<"\n";
        cout<<"  ops: ins="<<n_ins<<" del="<<n_del<<"\n";
        cout<<"  rebuilds="<<n_rebuild<<" fast="<<n_fast<<" slow="<<n_slow<<"\n";
    }
};

// ===== Dataset loading =====
struct TEdge { int u,v; long long ts; char op = '+'; };

vector<TEdge> load_dataset(const string& fn) {
    vector<TEdge> edges;
    ifstream fin(fn);
    if (!fin) { cerr<<"Cannot open "<<fn<<"\n"; return edges; }
    string line;
    while (getline(fin, line)) {
        if (line.empty()||line[0]=='#'||line[0]=='%') continue;
        istringstream iss(line);
        vector<string> tok;
        string s;
        while (iss >> s) tok.push_back(s);
        if (tok.size() < 3) continue;

        char op = '+';
        int u, v;
        long long ts;
        if (tok[0] == "+" || tok[0] == "I" || tok[0] == "i") {
            if (tok.size() < 4) continue;
            u = stoi(tok[1]); v = stoi(tok[2]); ts = stoll(tok[3]); op = '+';
        } else if (tok[0] == "-" || tok[0] == "D" || tok[0] == "d") {
            if (tok.size() < 4) continue;
            u = stoi(tok[1]); v = stoi(tok[2]); ts = stoll(tok[3]); op = '-';
        } else {
            u = stoi(tok[0]); v = stoi(tok[1]); ts = stoll(tok[2]);
            if (tok.size() >= 4 && (tok[3] == "-" || tok[3] == "D" || tok[3] == "d")) op = '-';
        }
        if (u != v) edges.push_back({u, v, ts, op});
    }
    sort(edges.begin(), edges.end(), [](auto& a, auto& b){ return a.ts<b.ts; });
    return edges;
}

void remap(vector<TEdge>& edges, int& nv) {
    unordered_map<int,int> m;
    int nxt = 0;
    for (auto& e : edges) {
        if (!m.count(e.u)) m[e.u]=nxt++;
        if (!m.count(e.v)) m[e.v]=nxt++;
        e.u=m[e.u]; e.v=m[e.v];
    }
    nv = nxt;
}

bool has_delete_events(const vector<TEdge>& edges) {
    for (const auto& e : edges) if (e.op == '-') return true;
    return false;
}

vector<pair<int,int>> active_edges_after(const vector<TEdge>& edges, int count) {
    EdgeSet active;
    for (int i = 0; i < count && i < (int)edges.size(); i++) {
        Edge e(edges[i].u, edges[i].v);
        if (edges[i].op == '-') active.erase(e);
        else active.insert(e);
    }
    vector<pair<int,int>> out;
    out.reserve(active.size());
    for (const auto& e : active) out.push_back({e.u, e.v});
    return out;
}

int main(int argc, char* argv[]) {
    if (argc<2) { cerr<<"Usage: "<<argv[0]<<" <dataset> [--validate]\n"; return 1; }
    string fn = argv[1];
    bool doVal = (argc>2 && string(argv[2])=="--validate");

    cout<<"========================================\n";
    cout<<"D-ESCENT: Dynamic EISTs\nDataset: "<<fn<<"\n";
    cout<<"========================================\n";

    Timer tm;
    tm.start();
    auto edges = load_dataset(fn);
    int nv = 0;
    remap(edges, nv);
    double loadT = tm.ms();
    cout<<"Loaded "<<edges.size()<<" temporal edges, "<<nv<<" vertices\n";
    cout<<"Load: "<<loadT<<" ms\n";

    bool fully_dynamic = has_delete_events(edges);
    vector<TEdge> ue;
    if (fully_dynamic) {
        ue = edges;
        cout<<"Dynamic events: "<<ue.size()<<" (insertions + deletions)\n";
    } else {
        EdgeSet seen;
        for (auto& e : edges) {
            Edge c(e.u,e.v);
            if (seen.insert(c).second) ue.push_back(e);
        }
        cout<<"Unique edges: "<<ue.size()<<"\n";
    }

    int initN = (int)(ue.size()*0.5);
    int insN  = (int)(ue.size()*0.3);
    int delN  = (int)(initN*0.2);

    // Phase 1: Init
    cout<<"\n--- Phase 1: Init ("<<initN<<" edges) ---\n";
    D_ESCENT ds;
    ds.init(nv);
    vector<pair<int,int>> initE;
    if (fully_dynamic) initE = active_edges_after(ue, initN);
    else for (int i=0;i<initN;i++) initE.push_back({ue[i].u, ue[i].v});
    tm.start();
    ds.build_initial(initE);
    double initT = tm.ms();
    cout<<"Init: "<<initT<<" ms\n";
    ds.print_stats();
    if (doVal) cout<<"Validation: "<<(ds.validate()?"PASS":"FAIL")<<"\n";

    // Phase 2: Insert
    cout<<"\n--- Phase 2: Insertions ("<<insN<<") ---\n";
    tm.start();
    for (int i=initN; i<initN+insN && i<(int)ue.size(); i++) {
        if (fully_dynamic && ue[i].op == '-') ds.delete_edge(ue[i].u, ue[i].v);
        else ds.insert_edge(ue[i].u, ue[i].v);
    }
    double insT = tm.ms();
    cout<<"Insert: "<<insT<<" ms (avg "<<(insN>0?insT/insN:0)<<" ms/op)\n";
    ds.print_stats();
    if (doVal) cout<<"Validation: "<<(ds.validate()?"PASS":"FAIL")<<"\n";

    // Phase 3: Delete
    cout<<"\n--- Phase 3: Deletions ("<<delN<<") ---\n";
    tm.start();
    if (fully_dynamic) {
        for (int i=initN+insN; i<initN+insN+delN && i<(int)ue.size(); i++) {
            if (ue[i].op == '-') ds.delete_edge(ue[i].u, ue[i].v);
            else ds.insert_edge(ue[i].u, ue[i].v);
        }
    } else {
        for (int i=initN-1; i>=initN-delN && i>=0; i--)
            ds.delete_edge(ue[i].u, ue[i].v);
    }
    double delT = tm.ms();
    cout<<"Delete: "<<delT<<" ms (avg "<<(delN>0?delT/delN:0)<<" ms/op)\n";
    ds.print_stats();
    if (doVal) cout<<"Validation: "<<(ds.validate()?"PASS":"FAIL")<<"\n";

    // Phase 4: Mixed
    int mixStart = fully_dynamic ? initN + insN + delN : initN + insN;
    int mixIns = (int)ue.size()-mixStart;
    int mixDel = fully_dynamic ? 0 : min(mixIns/3, delN);
    int mixTot = mixIns+mixDel;
    cout<<"\n--- Phase 4: Mixed ("<<mixIns<<" ins + "<<mixDel<<" del) ---\n";
    tm.start();
    int di=0;
    for (int i=mixStart; i<(int)ue.size(); i++) {
        if (fully_dynamic && ue[i].op == '-') ds.delete_edge(ue[i].u, ue[i].v);
        else ds.insert_edge(ue[i].u, ue[i].v);
        if (!fully_dynamic && di<mixDel && (i-initN-insN)%3==2) {
            int d=initN-delN-1-di;
            if (d>=0) { ds.delete_edge(ue[d].u, ue[d].v); di++; }
        }
    }
    double mixT = tm.ms();
    cout<<"Mixed: "<<mixT<<" ms (avg "<<(mixTot>0?mixT/mixTot:0)<<" ms/op)\n";
    ds.print_stats();
    if (doVal) cout<<"Validation: "<<(ds.validate()?"PASS":"FAIL")<<"\n";

    // Summary
    double dynT = insT+delT+mixT;
    int totOps = insN+delN+mixIns+mixDel;
    cout<<"\n========== SUMMARY ==========\n";
    cout<<"Dataset:       "<<fn<<"\n";
    cout<<"Total time:    "<<(initT+dynT)<<" ms\n";
    cout<<"  Init:        "<<initT<<" ms\n";
    cout<<"  Insertions:  "<<insT<<" ms ("<<insN<<" ops)\n";
    cout<<"  Deletions:   "<<delT<<" ms ("<<delN<<" ops)\n";
    cout<<"  Mixed:       "<<mixT<<" ms ("<<mixTot<<" ops)\n";
    cout<<"Dynamic ops:   "<<totOps<<"\n";
    cout<<"Avg time/op:   "<<(totOps>0?dynT/totOps:0)<<" ms\n";
    cout<<"Rebuilds:      "<<ds.n_rebuild<<"\n";
    cout<<"Fast path:     "<<ds.n_fast<<"\n";
    cout<<"Slow path:     "<<ds.n_slow<<"\n";
    cout<<"Final |V(H)|:  "<<ds.H.num_verts()<<"\n";
    cout<<"Final |E(H)|:  "<<ds.H.num_edges()<<"\n";
    cout<<"=============================\n";
    return 0;
}
