/*
 * Benchmark: Static vs Dynamic EIST Construction
 *
 * Five experiments measuring D-ESCENT dynamic EIST maintenance:
 *   Exp I   — Speedup vs batch size B ∈ {1, 10, 100, 1000, 10000}
 *   Exp II  — Fast-path fraction per round (B=100)
 *   Exp III — Tree quality: dynamic vs fresh static cost
 *   Exp IV  — Scaling: per-round cost as graph grows
 *   Exp V   — Rebuild trigger breakdown
 *
 * Uses sampled static checkpoints (MAX_SAMPLES) for tractability
 * when T rounds is large (B=1 or B=10).
 *
 * Compile: g++ -O3 -fopenmp -std=c++17 -o benchmark src/benchmark.cpp
 * Usage:   ./benchmark <dataset_file>
 */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <filesystem>
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
vector<int> largest_bicc_verts(const vector<vector<int>>& adj, const vector<bool>& active, int N) {
    vector<int> disc(N,-1), low(N,-1), par(N,-1);
    int timer_v = 0;
    vector<vector<int>> comps;
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
                if (low[v] >= disc[u]) pop_comp(u, v);
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
            if (!pushed && !fr.returning) stk.pop();
        }
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

// ===== Even-Tarjan st-numbering (OGDF pathfinder) =====
vector<int> compute_st_numbering(const Graph& H, int s, int t) {
    int N = H.N;
    auto verts = H.verts();
    int n = (int)verts.size();
    vector<int> g(N, 0);

    if (n <= 0) return g;
    if (n == 1) { g[s] = 1; return g; }
    if (n == 2) { g[s] = 1; g[t] = 2; return g; }

    vector<int> dfn(N, 0), low(N, 0), dfsParent(N, -1), followLow(N, -1);
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

    vector<bool> markedNode(N, false);
    markedNode[s] = markedNode[t] = true;
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
    vector<int> par1, par2;
    EdgeSet edg1, edg2;
    int root, t_vert;
};

TreePair build_trees(const Graph& H, const vector<int>& g, int root, int t_vert) {
    int N = H.N;
    TreePair tp;
    tp.root = root; tp.t_vert = t_vert;
    tp.par1.assign(N, -1); tp.par2.assign(N, -1);
    auto verts = H.verts();

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

    for (int v : verts) {
        if (v == root) continue;
        if (v == t_vert) {
            tp.par2[v] = root; tp.edg2.insert(Edge(v, root)); continue;
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
                dist[v] = dist[u]+1; q.push(v);
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
        if (H.adj[v1].empty()) continue;
        int tc = H.adj[v1][0];
        auto gc = compute_st_numbering(H, v1, tc);

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
    vector<int> bicc_lbl;
    int next_lbl = 0, main_lbl = -1;
    int TAU = 50;
    int MERGE_TAU = 512;
    double alpha = 1.5, growth_thr = 2.0;
    int n_ins=0, n_del=0, n_rebuild=0, n_fast=0, n_slow=0;
    EdgeSet non_tree;
    vector<EdgeSet> incident_non_tree;
    double cost_last = 0;
    int upd_ctr=0, merge_ctr=0, size_last=0;
    int g_size_last=0;
    bool pending_rebuild=false;

    // Per-round counters (reset each round for Exp II/V)
    int rnd_fast=0, rnd_slow=0;
    // Cumulative trigger counters (for Exp V)
    int trig_merge=0, trig_periodic=0, trig_pending=0;

    struct Replacement {
        Edge edge;
        int child = -1;
        int parent = -1;
        bool valid() const { return edge.valid() && child >= 0 && parent >= 0; }
    };

    void reset_round_counters() { rnd_fast = rnd_slow = 0; }

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
        if ((int)verts.size()<2) return;
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
        for (auto& e : H.edges)
            if (!tp.edg1.count(e) && !tp.edg2.count(e))
                add_non_tree_edge(e);

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
            n_fast++; rnd_fast++;
        } else {
            merge_ctr++;
            bool lazy = (merge_ctr<MERGE_TAU) && (G.num_verts()<(int)(growth_thr*g_size_last));
            if (lazy) { pending_rebuild=true; n_fast++; rnd_fast++; }
            else { do_full_rebuild(5); n_slow++; rnd_slow++; trig_merge++; }
        }
        periodic();
    }

    void delete_edge(int u, int v) {
        n_del++;
        if (!G.has_edge(u,v)) return;
        G.remove_edge(u,v); upd_ctr++;
        Edge e(u,v);

        if (!H.has_edge(u,v)) { n_fast++; rnd_fast++; periodic(); return; }
        H.remove_edge(u,v);
        bool i1 = tp.edg1.count(e), i2 = tp.edg2.count(e);

        if (!i1 && !i2) {
            remove_non_tree_edge(e); n_fast++; rnd_fast++;
        } else if (i1 && !i2) {
            tp.edg1.erase(e);
            int child = -1;
            if (tp.par1[u]==v) child=u;
            else if (tp.par1[v]==u) child=v;
            else child = (g[u]>g[v])?u:v;
            Replacement r = find_repl(child, true);
            if (r.valid()) {
                tp.par1[r.child] = r.parent;
                tp.edg1.insert(r.edge); remove_non_tree_edge(r.edge); n_fast++; rnd_fast++;
            } else { rebuild_H(1,root); n_slow++; rnd_slow++; }
        } else if (!i1 && i2) {
            tp.edg2.erase(e);
            int child = -1;
            if (tp.par2[u]==v) child=u;
            else if (tp.par2[v]==u) child=v;
            else child = (g[u]<g[v])?u:v;
            Replacement r = find_repl(child, false);
            if (r.valid()) {
                tp.par2[r.child] = r.parent;
                tp.edg2.insert(r.edge); remove_non_tree_edge(r.edge); n_fast++; rnd_fast++;
            } else { rebuild_H(1,root); n_slow++; rnd_slow++; }
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
                n_fast++; rnd_fast++;
            } else { rebuild_H(1,root); n_slow++; rnd_slow++; }
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
                trig_pending++;
            }
            return;
        }
        auto verts = H.verts();
        if ((int)verts.size()>=2) {
            double cc = compute_cost(tp, verts, H.N);
            if (cc > alpha*cost_last) { rebuild_H(5); trig_periodic++; }
            else cost_last = cc;
        }
    }

    double current_cost() {
        auto verts = H.verts();
        if ((int)verts.size()<2) return 0;
        return compute_cost(tp, verts, H.N);
    }
};

// ===== Static EIST reconstruction with phasewise timing =====
struct StaticResult {
    double total_ms=0, bicc_ms=0, root_ms=0, st_ms=0, tree_ms=0;
    double cost=0;
    int bicc_size=0;
};

// ===== Dataset loading =====
struct TEdge { int u,v; long long ts; char op = '+'; };

// use_count: how many events from the start of all_events to use
StaticResult static_rebuild_full(const vector<TEdge>& all_events,
                                  int use_count, int nv, int K=32) {
    StaticResult r;
    Timer tm;

    tm.start();
    Graph G; G.init(nv);
    for (int i = 0; i < use_count && i < (int)all_events.size(); i++) {
        const auto& ev = all_events[i];
        if (ev.op == '-') G.remove_edge(ev.u, ev.v);
        else G.add_edge(ev.u, ev.v);
    }
    Graph H = extract_largest_bicc(G);
    auto verts = H.verts();
    r.bicc_ms = tm.ms();
    r.bicc_size = (int)verts.size();
    if ((int)verts.size() < 2) { r.total_ms = r.bicc_ms; return r; }

    tm.start();
    int root = double_sweep_bfs(H, K);
    if (root < 0 || root >= H.N || H.adj[root].empty()) { r.total_ms = r.bicc_ms; return r; }
    int t_vert = H.adj[root][0];
    r.root_ms = tm.ms();

    tm.start();
    auto g = compute_st_numbering(H, root, t_vert);
    r.st_ms = tm.ms();

    tm.start();
    auto tp = build_trees(H, g, root, t_vert);
    r.tree_ms = tm.ms();

    r.total_ms = r.bicc_ms + r.root_ms + r.st_ms + r.tree_ms;
    r.cost = compute_cost(tp, verts, H.N);
    return r;
}

// ===== Sampling utility =====
vector<int> sample_indices(int T, int max_samples) {
    vector<int> idx;
    if (T <= max_samples) {
        idx.resize(T);
        iota(idx.begin(), idx.end(), 0);
        return idx;
    }
    set<int> s;
    for (int i = 0; i < max_samples; i++) {
        int j = (int)((long long)(T-1) * i / (max_samples-1));
        s.insert(j);
    }
    idx.assign(s.begin(), s.end());
    return idx;
}

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

// ===== Benchmark constants =====
const int B_VALUES[] = {1, 10, 100, 1000, 10000};
const int NUM_B = 5;
const double INIT_FRAC = 0.5;
const int K_STATIC = 32;
const int MAX_SAMPLES = 300;
const int MAX_POOL = 50000;         // cap pool to limit runtime on large datasets
const int B_MAIN = 100;             // primary batch size for single-pass recording
const int TIMEOUT_PER_DS = 1800;    // 30 min per dataset

// ===== Per-round data collected during the single dynamic pass =====
struct RoundData {
    double dyn_ms;         // dynamic time for this batch
    int rnd_fast, rnd_slow;
    double dyn_cost;       // cost of dynamic trees after this batch
    int n_edges;           // total accumulated edges after this batch
    int cum_rebuilds;
    int cum_trig_merge, cum_trig_periodic, cum_trig_pending;
};

// ===== Main benchmark =====
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dataset_file> [--out-dir DIR]" << endl;
        return 1;
    }

    string fn = argv[1];
    string out_dir = "results";
    for (int i = 2; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--out-dir" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            cerr << "Unknown argument: " << arg << endl;
            return 1;
        }
    }

    cout << "========================================" << endl;
    cout << "D-ESCENT Benchmark Suite (single-pass)" << endl;
    cout << "Dataset: " << fn << endl;
    cout << "Output: " << out_dir << endl;
#ifdef _OPENMP
    cout << "OpenMP threads: " << omp_get_max_threads() << endl;
#else
    cout << "OpenMP: disabled" << endl;
#endif
    cout << "========================================" << endl;

    // Load and remap dataset
    auto edges = load_dataset(fn);
    int nv = 0;
    remap(edges, nv);
    cout << "Loaded " << edges.size() << " temporal edges, " << nv << " vertices" << endl;

    bool fully_dynamic = has_delete_events(edges);
    vector<TEdge> ue;
    if (fully_dynamic) {
        ue = edges;
        cout << "Dynamic events: " << ue.size() << " (insertions + deletions)" << endl;
    } else {
        EdgeSet seen;
        for (auto& e : edges) {
            Edge c(e.u, e.v);
            if (seen.insert(c).second) ue.push_back(e);
        }
        cout << "Unique edges: " << ue.size() << endl;
    }

    // Split into init and update pool. For fully dynamic streams, init is
    // the active graph after replaying the first INIT_FRAC events.
    int init_sz = (int)(ue.size() * INIT_FRAC);
    vector<pair<int,int>> init_edges;
    vector<TEdge> pool;
    if (fully_dynamic) {
        init_edges = active_edges_after(ue, init_sz);
    } else {
        for (int i = 0; i < init_sz; i++)
            init_edges.push_back({ue[i].u, ue[i].v});
    }
    for (int i = init_sz; i < (int)ue.size(); i++)
        pool.push_back(ue[i]);

    // Cap pool for large datasets
    if ((int)pool.size() > MAX_POOL) {
        cout << "Pool capped: " << pool.size() << " -> " << MAX_POOL << endl;
        pool.resize(MAX_POOL);
    }

    cout << "Init edges: " << init_edges.size()
         << ", Pool edges: " << pool.size() << endl;

    // Extract dataset name
    string ds_name = fn;
    size_t last_slash = ds_name.find_last_of("/\\");
    if (last_slash != string::npos) ds_name = ds_name.substr(last_slash + 1);
    size_t dot = ds_name.find_last_of('.');
    if (dot != string::npos) ds_name = ds_name.substr(0, dot);

    filesystem::create_directories(out_dir);

    // Build event list for static rebuilds.
    vector<TEdge> accum;
    if (fully_dynamic) {
        accum = ue;
    } else {
        accum.reserve(ue.size());
        for (const auto& e : ue) accum.push_back(e);
    }

    int pool_sz = (int)pool.size();
    int T = pool_sz / B_MAIN;
    if (T < 1) { cerr << "Pool too small for B=" << B_MAIN << endl; return 1; }

    // =========================================================
    // PHASE 1: Single dynamic pass (B=100 granularity)
    // =========================================================
    cout << "\n--- Dynamic Pass (" << pool_sz << " edges, "
         << T << " rounds of " << B_MAIN << ") ---" << endl;

    D_ESCENT ds;
    ds.init(nv);
    ds.build_initial(init_edges);

    vector<RoundData> rounds(T);
    Timer total_wall; total_wall.start();
    Timer dyn_wall; dyn_wall.start();

    for (int r = 0; r < T; r++) {
        int start = r * B_MAIN;
        int end = min(start + B_MAIN, pool_sz);

        ds.reset_round_counters();
        Timer dt; dt.start();
        for (int i = start; i < end; i++) {
            if (pool[i].op == '-') ds.delete_edge(pool[i].u, pool[i].v);
            else ds.insert_edge(pool[i].u, pool[i].v);
        }

        rounds[r].dyn_ms = dt.ms();
        rounds[r].rnd_fast = ds.rnd_fast;
        rounds[r].rnd_slow = ds.rnd_slow;
        rounds[r].dyn_cost = ds.current_cost();
        rounds[r].n_edges = ds.G.num_edges();
        rounds[r].cum_rebuilds = ds.n_rebuild;
        rounds[r].cum_trig_merge = ds.trig_merge;
        rounds[r].cum_trig_periodic = ds.trig_periodic;
        rounds[r].cum_trig_pending = ds.trig_pending;

        if (r % 100 == 0)
            cout << "  round " << r << "/" << T
                 << " | edges=" << rounds[r].n_edges
                 << " | dyn=" << fixed << setprecision(2) << rounds[r].dyn_ms << "ms"
                 << " | rebuilds=" << ds.n_rebuild << endl;
    }

    double dyn_total_ms = dyn_wall.ms();
    double fast_pct = (ds.n_fast + ds.n_slow > 0) ?
        100.0 * ds.n_fast / (ds.n_fast + ds.n_slow) : 0;

    cout << fixed << setprecision(1)
         << "Dynamic total: " << dyn_total_ms/1000 << "s"
         << " | rebuilds=" << ds.n_rebuild
         << " | fast=" << fast_pct << "%" << endl;

    // =========================================================
    // PHASE 2: Static sampling at checkpoints
    // =========================================================
    auto samp_idx = sample_indices(T, MAX_SAMPLES);
    int S = (int)samp_idx.size();
    unordered_set<int> samp_set(samp_idx.begin(), samp_idx.end());

    cout << "\n--- Static Sampling (" << S << " checkpoints) ---" << endl;

    struct StaticSample {
        int round;
        StaticResult sr;
    };
    vector<StaticSample> statics(S);
    Timer sta_wall; sta_wall.start();

    for (int si = 0; si < S; si++) {
        int r = samp_idx[si];
        int ne = init_sz + min((r + 1) * B_MAIN, pool_sz);
        statics[si].round = r;
        statics[si].sr = static_rebuild_full(accum, ne, nv, K_STATIC);

        if (si % 50 == 0)
            cout << "  sample " << si << "/" << S
                 << " | edges=" << ne
                 << " | static=" << fixed << setprecision(2)
                 << statics[si].sr.total_ms << "ms" << endl;
    }
    cout << "Static sampling: " << fixed << setprecision(1)
         << sta_wall.ms()/1000 << "s" << endl;

    // =========================================================
    // PHASE 3: Write all 5 experiment CSVs
    // =========================================================
    cout << "\n--- Writing Results ---" << endl;

    // ----- Exp I: Speedup vs Batch Size -----
    {
        string csv = (filesystem::path(out_dir) / ("exp1_speedup_" + ds_name + ".csv")).string();
        ofstream out(csv);
        out << "B,T_rounds,S_samples,dyn_total_ms,sta_total_ms,speedup,"
            << "sta_bicc_avg_ms,sta_root_avg_ms,sta_st_avg_ms,sta_tree_avg_ms,"
            << "dyn_rebuilds,dyn_fast_pct" << endl;

        for (int bi = 0; bi < NUM_B; bi++) {
            int B = B_VALUES[bi];
            int T_B = pool_sz / B;
            if (T_B < 1) continue;

            // For each B, compute sampled static cost at the B's round granularity
            // Round r at batch size B has init_sz + (r+1)*B accumulated edges
            auto samps_B = sample_indices(T_B, MAX_SAMPLES);
            int S_B = (int)samps_B.size();

            vector<double> sampled_sta_times(S_B, 0);
            double bicc_sum=0, root_sum=0, st_sum=0, tree_sum=0;

            for (int si = 0; si < S_B; si++) {
                int r_B = samps_B[si];
                int ne = init_sz + (r_B + 1) * B;
                if (ne > (int)accum.size()) ne = (int)accum.size();
                // Find closest pre-computed static sample or compute fresh
                // Check if we already have this edge count from the B=100 sampling
                bool found = false;
                for (int j = 0; j < S; j++) {
                    int sampled_events = init_sz + min((statics[j].round + 1) * B_MAIN, pool_sz);
                    if (sampled_events == ne) {
                        sampled_sta_times[si] = statics[j].sr.total_ms;
                        bicc_sum += statics[j].sr.bicc_ms;
                        root_sum += statics[j].sr.root_ms;
                        st_sum += statics[j].sr.st_ms;
                        tree_sum += statics[j].sr.tree_ms;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    // Interpolate from nearest pre-computed samples
                    // Find bracketing static samples by edge count
                    int best_lo = -1, best_hi = -1;
                    for (int j = 0; j < S; j++) {
                        int sne = init_sz + min((statics[j].round + 1) * B_MAIN, pool_sz);
                        int best_lo_events = best_lo < 0 ? -1 : init_sz + min((statics[best_lo].round + 1) * B_MAIN, pool_sz);
                        int best_hi_events = best_hi < 0 ? INT_MAX : init_sz + min((statics[best_hi].round + 1) * B_MAIN, pool_sz);
                        if (sne <= ne && (best_lo < 0 || sne > best_lo_events))
                            best_lo = j;
                        if (sne >= ne && (best_hi < 0 || sne < best_hi_events))
                            best_hi = j;
                    }
                    if (best_lo >= 0 && best_hi >= 0 && best_lo != best_hi) {
                        int lo_ne = init_sz + min((statics[best_lo].round + 1) * B_MAIN, pool_sz);
                        int hi_ne = init_sz + min((statics[best_hi].round + 1) * B_MAIN, pool_sz);
                        double frac = (double)(ne - lo_ne) / max(1, hi_ne - lo_ne);
                        sampled_sta_times[si] = statics[best_lo].sr.total_ms +
                            frac * (statics[best_hi].sr.total_ms - statics[best_lo].sr.total_ms);
                        double lo_bicc = statics[best_lo].sr.bicc_ms, hi_bicc = statics[best_hi].sr.bicc_ms;
                        bicc_sum += lo_bicc + frac * (hi_bicc - lo_bicc);
                        root_sum += statics[best_lo].sr.root_ms + frac * (statics[best_hi].sr.root_ms - statics[best_lo].sr.root_ms);
                        st_sum += statics[best_lo].sr.st_ms + frac * (statics[best_hi].sr.st_ms - statics[best_lo].sr.st_ms);
                        tree_sum += statics[best_lo].sr.tree_ms + frac * (statics[best_hi].sr.tree_ms - statics[best_lo].sr.tree_ms);
                    } else if (best_lo >= 0) {
                        sampled_sta_times[si] = statics[best_lo].sr.total_ms;
                        bicc_sum += statics[best_lo].sr.bicc_ms;
                        root_sum += statics[best_lo].sr.root_ms;
                        st_sum += statics[best_lo].sr.st_ms;
                        tree_sum += statics[best_lo].sr.tree_ms;
                    } else if (best_hi >= 0) {
                        sampled_sta_times[si] = statics[best_hi].sr.total_ms;
                        bicc_sum += statics[best_hi].sr.bicc_ms;
                        root_sum += statics[best_hi].sr.root_ms;
                        st_sum += statics[best_hi].sr.st_ms;
                        tree_sum += statics[best_hi].sr.tree_ms;
                    }
                }
            }

            // Interpolate total static via trapezoidal rule
            double sta_total = 0;
            if (S_B == T_B) {
                for (int i = 0; i < S_B; i++) sta_total += sampled_sta_times[i];
            } else if (S_B >= 2) {
                for (int i = 0; i < S_B - 1; i++) {
                    int gap = samps_B[i+1] - samps_B[i];
                    sta_total += gap * (sampled_sta_times[i] + sampled_sta_times[i+1]) / 2.0;
                }
            } else if (S_B == 1) {
                sta_total = sampled_sta_times[0] * T_B;
            }

            double speedup = (dyn_total_ms > 0) ? sta_total / dyn_total_ms : 0;

            out << B << "," << T_B << "," << S_B << ","
                << fixed << setprecision(2) << dyn_total_ms << ","
                << sta_total << "," << setprecision(2) << speedup << ","
                << (S_B>0 ? bicc_sum/S_B : 0) << ","
                << (S_B>0 ? root_sum/S_B : 0) << ","
                << (S_B>0 ? st_sum/S_B : 0) << ","
                << (S_B>0 ? tree_sum/S_B : 0) << ","
                << ds.n_rebuild << ","
                << setprecision(1) << fast_pct << endl;

            cout << "  Exp I: B=" << B << " T=" << T_B
                 << " speedup=" << setprecision(2) << speedup << "x" << endl;
        }
        out.close();
        cout << "  -> " << csv << endl;
    }

    // ----- Exp II: Fast-path Breakdown -----
    {
        string csv = (filesystem::path(out_dir) / ("exp2_fastpath_" + ds_name + ".csv")).string();
        ofstream out(csv);
        out << "round,n_edges,fast,slow,fast_pct" << endl;
        for (int r = 0; r < T; r++) {
            int total_ops = rounds[r].rnd_fast + rounds[r].rnd_slow;
            double fpct = total_ops > 0 ? 100.0*rounds[r].rnd_fast/total_ops : 100.0;
            out << r << "," << rounds[r].n_edges << ","
                << rounds[r].rnd_fast << "," << rounds[r].rnd_slow << ","
                << fixed << setprecision(1) << fpct << endl;
        }
        out.close();
        cout << "  -> " << csv << " (Exp II, " << T << " rounds)" << endl;
    }

    // ----- Exp III: Tree Quality -----
    {
        string csv = (filesystem::path(out_dir) / ("exp3_quality_" + ds_name + ".csv")).string();
        ofstream out(csv);
        out << "round,n_edges,dyn_cost,sta_cost,ratio" << endl;
        for (int si = 0; si < S; si++) {
            int r = statics[si].round;
            double dyn_cost = rounds[r].dyn_cost;
            double sta_cost = statics[si].sr.cost;
            double ratio = (sta_cost > 0) ? dyn_cost / sta_cost : 1.0;
            out << r << "," << rounds[r].n_edges << ","
                << fixed << setprecision(2) << dyn_cost << ","
                << sta_cost << "," << setprecision(4) << ratio << endl;
        }
        out.close();
        cout << "  -> " << csv << " (Exp III, " << S << " samples)" << endl;
    }

    // ----- Exp IV: Scaling with Graph Size -----
    {
        string csv = (filesystem::path(out_dir) / ("exp4_scaling_" + ds_name + ".csv")).string();
        ofstream out(csv);
        out << "round,n_edges,bicc_size,dyn_ms,sta_ms,"
            << "sta_bicc_ms,sta_root_ms,sta_st_ms,sta_tree_ms" << endl;
        for (int si = 0; si < S; si++) {
            int r = statics[si].round;
            auto& sr = statics[si].sr;
            out << r << "," << rounds[r].n_edges << ","
                << sr.bicc_size << ","
                << fixed << setprecision(3) << rounds[r].dyn_ms << ","
                << sr.total_ms << ","
                << sr.bicc_ms << "," << sr.root_ms << ","
                << sr.st_ms << "," << sr.tree_ms << endl;
        }
        out.close();
        cout << "  -> " << csv << " (Exp IV, " << S << " samples)" << endl;
    }

    // ----- Exp V: Rebuild Triggers -----
    {
        string csv = (filesystem::path(out_dir) / ("exp5_triggers_" + ds_name + ".csv")).string();
        ofstream out(csv);
        out << "round,n_edges,cum_rebuilds,cum_trig_merge,cum_trig_periodic,cum_trig_pending" << endl;
        int report_every = max(1, T / 300);
        for (int r = 0; r < T; r++) {
            if (r % report_every == 0 || r == T-1) {
                out << r << "," << rounds[r].n_edges << ","
                    << rounds[r].cum_rebuilds << ","
                    << rounds[r].cum_trig_merge << ","
                    << rounds[r].cum_trig_periodic << ","
                    << rounds[r].cum_trig_pending << endl;
            }
        }
        out.close();
        cout << "  -> " << csv << " (Exp V, " << T << " rounds)" << endl;
    }

    // =========================================================
    // Summary
    // =========================================================
    cout << "\n========== COMPLETE ==========" << endl;
    cout << fixed << setprecision(1);
    cout << "Total wall time: " << total_wall.ms()/1000.0 << " s" << endl;
    cout << "Dynamic time:    " << dyn_total_ms/1000.0 << " s" << endl;
    cout << "Static sampling: " << sta_wall.ms()/1000.0 << " s" << endl;
    cout << "Rebuilds: " << ds.n_rebuild
         << " (merge=" << ds.trig_merge
         << ", periodic=" << ds.trig_periodic
         << ", pending=" << ds.trig_pending << ")" << endl;
    cout << "Fast path: " << setprecision(1) << fast_pct << "%" << endl;
    cout << "==============================" << endl;

    return 0;
}
