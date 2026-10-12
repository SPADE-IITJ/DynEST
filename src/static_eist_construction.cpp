#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <stack>
#include <queue>
#include <list>
#include <map>


using namespace std;

struct Graph {
    int n;
    vector<vector<int>> adj;
    
    Graph(int vertices) : n(vertices), adj(vertices) {}
    
    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

class BiCCExtractor {
private:
    const Graph& g;
    vector<int> disc, low, parent;
    vector<bool> visited;
    int timer;

    stack<pair<int,int>> edge_stack;

public:
    BiCCExtractor(const Graph& graph)
        : g(graph),
          disc(graph.n, -1),
          low(graph.n, -1),
          parent(graph.n, -1),
          visited(graph.n, false),
          timer(0) {}

    Graph extract_largest() {
        vector<pair<int,int>> best_bicc;
        int best_size = 0;

        for (int start = 0; start < g.n; start++) {
            if (visited[start]) continue;

            stack<int> dfs_stack;
            stack<int> iter_stack;

            dfs_stack.push(start);
            parent[start] = -1;

            while (!dfs_stack.empty()) {
                int u = dfs_stack.top();

                if (!visited[u]) {
                    visited[u] = true;
                    disc[u] = low[u] = timer++;
                    iter_stack.push(0);
                }

                bool finished = true;

                while (!iter_stack.empty()) {
                    int& idx = iter_stack.top();
                    if (idx < (int)g.adj[u].size()) {
                        int v = g.adj[u][idx++];
                        if (!visited[v]) {
                            parent[v] = u;
                            edge_stack.push({u, v});
                            dfs_stack.push(v);
                            finished = false;
                            break;
                        } else if (v != parent[u] && disc[v] < disc[u]) {
                            edge_stack.push({u, v});
                            low[u] = min(low[u], disc[v]);
                        }
                    } else {
                        iter_stack.pop();
                        break;
                    }
                }

                if (!finished) continue;

                dfs_stack.pop();

                int p = parent[u];
                if (p != -1) {
                    low[p] = min(low[p], low[u]);

                    if (low[u] >= disc[p]) {
                        unordered_set<int> verts;
                        vector<pair<int,int>> curr;

                        while (!edge_stack.empty()) {
                            auto e = edge_stack.top();
                            edge_stack.pop();
                            curr.push_back(e);
                            verts.insert(e.first);
                            verts.insert(e.second);
                            if (e.first == p && e.second == u) break;
                        }

                        if ((int)verts.size() > best_size) {
                            best_size = verts.size();
                            best_bicc.swap(curr);
                        }
                    }
                }
            }

            // leftover edges
            if (!edge_stack.empty()) {
                unordered_set<int> verts;
                vector<pair<int,int>> curr;
                while (!edge_stack.empty()) {
                    auto e = edge_stack.top();
                    edge_stack.pop();
                    curr.push_back(e);
                    verts.insert(e.first);
                    verts.insert(e.second);
                }
                if ((int)verts.size() > best_size) {
                    best_size = verts.size();
                    best_bicc.swap(curr);
                }
            }
        }

        if (best_bicc.empty()) return Graph(0);

        // Build graph from best BiCC
        unordered_map<int,int> remap;
        int id = 0;
        for (auto& e : best_bicc) {
            if (!remap.count(e.first)) remap[e.first] = id++;
            if (!remap.count(e.second)) remap[e.second] = id++;
        }

        Graph h(remap.size());
        set<pair<int,int>> seen;

        for (auto& e : best_bicc) {
            int u = remap[e.first];
            int v = remap[e.second];
            auto key = minmax(u, v);
            if (seen.insert(key).second) {
                h.add_edge(u, v);
            }
        }

        return h;
    }
};


// ST-numbering using OGDF-style pathfinder (Even-Tarjan algorithm)
// Precondition: graph is biconnected, s and t are adjacent.
// Returns g[] with g[s]=1, g[t]=n.
vector<int> compute_st_numbering(const Graph& graph, int s, int t) {
    int n = graph.n;
    if (n < 2) return {};
    if (n == 2) { vector<int> r(2,0); r[s]=1; r[t]=2; return r; }

    vector<int> dfn(n, 0), low_v(n, 0), dfs_par(n, -1), follow_low(n, -1);
    int cnt = 1;
    dfn[t] = cnt++; low_v[t] = dfn[t];
    dfn[s] = cnt++; low_v[s] = dfn[s];

    {
        struct Fr { int u; int idx; };
        stack<Fr> stk;
        stk.push({s, 0});
        while (!stk.empty()) {
            auto& fr = stk.top();
            int u = fr.u;
            bool pushed = false;
            while (fr.idx < (int)graph.adj[u].size()) {
                int w = graph.adj[u][fr.idx]; fr.idx++;
                if (dfn[w] == 0) {
                    dfs_par[w] = u;
                    dfn[w] = cnt++; low_v[w] = dfn[w];
                    stk.push({w, 0});
                    pushed = true;
                    break;
                } else if (low_v[u] > dfn[w]) {
                    low_v[u] = dfn[w];
                    follow_low[u] = w;
                }
            }
            if (!pushed) {
                stk.pop();
                int p = dfs_par[u];
                if (p >= 0 && low_v[p] > low_v[u]) {
                    low_v[p] = low_v[u];
                    follow_low[p] = u;
                }
            }
        }
    }
    if (low_v[t] > low_v[s]) low_v[t] = low_v[s];

    // Pathfinder phase
    vector<bool> marked_n(n, false);
    marked_n[s] = marked_n[t] = true;
    set<pair<int,int>> marked_e;
    auto mke = [](int a, int b) { return make_pair(min(a,b), max(a,b)); };
    marked_e.insert(mke(s, t));

    vector<int> adj_idx(n, 0);
    stack<int> ns;
    ns.push(t);
    ns.push(s);

    vector<int> result(n, 0);
    int num = 1;
    int v = ns.top(); ns.pop();

    while (v != t) {
        bool found = false;
        while (adj_idx[v] < (int)graph.adj[v].size()) {
            int w = graph.adj[v][adj_idx[v]]; adj_idx[v]++;
            if (marked_e.count(mke(v, w))) continue;
            marked_e.insert(mke(v, w));

            if (dfs_par[w] == v) {
                vector<int> path;
                path.push_back(v);
                while (!marked_n[w]) {
                    path.push_back(w);
                    marked_n[w] = true;
                    int nxt = follow_low[w];
                    if (nxt < 0) break;
                    marked_e.insert(mke(w, nxt));
                    w = nxt;
                }
                for (int i = (int)path.size()-1; i >= 0; i--)
                    ns.push(path[i]);
                found = true;
                break;
            } else if (dfn[v] < dfn[w]) {
                vector<int> path;
                path.push_back(v);
                while (!marked_n[w]) {
                    path.push_back(w);
                    marked_n[w] = true;
                    int par = dfs_par[w];
                    if (par < 0) break;
                    marked_e.insert(mke(w, par));
                    w = par;
                }
                for (int i = (int)path.size()-1; i >= 0; i--)
                    ns.push(path[i]);
                found = true;
                break;
            }
        }
        if (!found) {
            result[v] = num++;
            adj_idx[v] = 0;
        }
        v = ns.top(); ns.pop();
    }
    result[t] = num++;
    return result;
}



// Itai-Rodeh EIST construction
class ItaiRodeh {
private:
    const Graph& g;
    int root;
    int t_vert;
    vector<int> st_num;
    
public:
    ItaiRodeh(const Graph& graph, int r, int sink, const vector<int>& numbering)
        : g(graph), root(r), t_vert(sink), st_num(numbering) {}
    
    pair<vector<int>, vector<int>> construct() {
        int n = g.n;
        int t = t_vert;
        vector<int> parS(n, -1);
        vector<int> parT(n, -1);
        
        // Tree S: low-parent tree
        for (int v = 0; v < n; v++) {
            if (v == root) continue;
            
            int best_parent = -1;
            for (int u : g.adj[v]) {
                if (st_num[u] < st_num[v]) {
                    if (v == t && u == root) continue;
                    if (best_parent == -1 || st_num[u] > st_num[best_parent]) {
                        best_parent = u;
                    }
                }
            }
            parS[v] = best_parent;
        }
        
        // Tree T: high-parent tree
        parT[t] = root;
        for (int v = 0; v < n; v++) {
            if (v == root || v == t) continue;
            
            int best_parent = -1;
            for (int u : g.adj[v]) {
                if (st_num[u] > st_num[v]) {
                    if (best_parent == -1 || st_num[u] < st_num[best_parent]) {
                        best_parent = u;
                    }
                }
            }
            parT[v] = best_parent;
        }
        
        return {parS, parT};
    }
};

class DoubleSweepBFS {
private:
    const Graph& g;
    
    // BFS returning farthest vertex and its distance
    pair<int, int> bfs_farthest(int start) {
        int n = g.n;
        vector<int> dist(n, -1);
        queue<int> q;
        
        q.push(start);
        dist[start] = 0;
        
        int farthest = start;
        int max_dist = 0;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (int v : g.adj[u]) {
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    q.push(v);
                    
                    if (dist[v] > max_dist) {
                        max_dist = dist[v];
                        farthest = v;
                    }
                }
            }
        }
        
        return {farthest, max_dist};
    }
    
    // BFS that records parent pointers for path reconstruction
    pair<int, vector<int>> bfs_farthest_with_path(int start) {
        int n = g.n;
        vector<int> dist(n, -1);
        vector<int> parent(n, -1);
        queue<int> q;
        
        q.push(start);
        dist[start] = 0;
        
        int farthest = start;
        int max_dist = 0;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (int v : g.adj[u]) {
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    parent[v] = u;
                    q.push(v);
                    
                    if (dist[v] > max_dist) {
                        max_dist = dist[v];
                        farthest = v;
                    }
                }
            }
        }
        
        return {farthest, parent};
    }
    
    // Reconstruct path from parent pointers
    vector<int> reconstruct_path(int v1, int v2, const vector<int>& parent) {
        vector<int> path;
        int curr = v2;
        
        while (curr != -1) {
            path.push_back(curr);
            if (curr == v1) break;
            curr = parent[curr];
        }
        
        reverse(path.begin(), path.end());
        return path;
    }
    
public:
    DoubleSweepBFS(const Graph& graph) : g(graph) {}
    
    // Single double-sweep: returns middle vertex/vertices of diameter path
    vector<int> single_sweep(int seed) {
        // First BFS: seed → v1
        auto [v1, dist1] = bfs_farthest(seed);
        
        // Second BFS: v1 → v2, with path reconstruction
        auto [v2, parent] = bfs_farthest_with_path(v1);
        
        // Reconstruct path P_{v1,v2}
        vector<int> path = reconstruct_path(v1, v2, parent);
        
        // Extract middle vertex/vertices
        vector<int> middle;
        int len = path.size();
        
        if (len == 0) return middle;
        
        if (len % 2 == 1) {
            // Odd length: single middle
            middle.push_back(path[len / 2]);
        } 
        else {
            // Even length: two middles
            middle.push_back(path[len / 2 - 1]);
            // middle.push_back(path[len / 2]);
        }
        
        return middle;
    }
};


// Depth metrics computation — O(N) BFS-based
pair<int, double> compute_depth_metrics(const vector<int>& parent, int root) {
    int n = parent.size();
    vector<vector<int>> children(n);
    for (int v = 0; v < n; v++) {
        if (v != root && parent[v] >= 0 && parent[v] < n)
            children[parent[v]].push_back(v);
    }
    vector<int> depth(n, 0);
    int max_depth = 0;
    double avg_depth = 0.0;
    int count = 0;
    queue<int> bfs;
    bfs.push(root);
    while (!bfs.empty()) {
        int u = bfs.front(); bfs.pop();
        for (int c : children[u]) {
            depth[c] = depth[u] + 1;
            if (depth[c] > max_depth) max_depth = depth[c];
            avg_depth += depth[c]; count++;
            bfs.push(c);
        }
    }
    
    if (count > 0) {
        avg_depth /= count;
    }
    
    return {max_depth, avg_depth};
}


unordered_set<int> generate_candidate_roots(const Graph& h, int K = 32) {
    DoubleSweepBFS sweeper(h);
    unordered_set<int> candidates;
    
    for(int i = 0; i < K; i++) {
        int start = i % h.n;

        vector<int> middles = sweeper.single_sweep(start);

        for (int v : middles) {
            candidates.insert(v);
        }

    }
    
    cout << "Generated " << candidates.size() 
         << " candidate roots from " << K << " sweeps" << endl;
    
    return candidates;
}

double compute_cost(const vector<int>& parS, const vector<int>& parT,
                   int root) {
    auto [S_max, S_avg] = compute_depth_metrics(parS, root);
    auto [T_max, T_avg] = compute_depth_metrics(parT, root);
    
    return max(S_max, T_max) + S_avg + T_avg;
}


struct EISTResult {
    vector<int> parS;
    vector<int> parT;
    int root;
    int neighbor;
    double cost;
    double cand_ms;
    double st_ms;
    double eist_ms;
};

EISTResult root_aware_eist(const Graph& h) {
    if (h.n == 0) {
        return {vector<int>(), vector<int>(), -1, -1, 0.0, 0.0, 0.0, 0.0};
    }
    
    if (h.n == 1) {
        return {vector<int>(1, -1), vector<int>(1, -1), 0, -1, 0.0, 0.0, 0.0, 0.0};
    }
    
    auto cs = chrono::high_resolution_clock::now();
    const int K = 32;
    unordered_set<int> candidates = generate_candidate_roots(h, K);
    candidates.insert(0);
    
    if (candidates.empty()) {
        candidates.insert(0);
    }
    auto ce = chrono::high_resolution_clock::now();
    double cand_ms = chrono::duration<double, milli>(ce - cs).count();
    
    double best_cost = numeric_limits<double>::max();
    vector<int> best_parS, best_parT;
    int best_root = -1;
    int best_neighbor = -1;
    double st_ms = 0.0;
    double eist_ms = 0.0;
    
    for (int r : candidates) {
        if (h.adj[r].empty()) {
            continue;
        }
        
        // Try all neighbors of root until one produces valid st-numbering
        int t = -1;
        vector<int> st_num;
        for (int candidate_t : h.adj[r]) {
            auto s1 = chrono::high_resolution_clock::now();
            st_num = compute_st_numbering(h, r, candidate_t);
            auto s2 = chrono::high_resolution_clock::now();
            st_ms += chrono::duration<double, milli>(s2 - s1).count();
            
            // Check if this neighbor works
            if (!st_num.empty() && st_num[r] == 1 && st_num[candidate_t] == h.n) {
                t = candidate_t;
                break;  // Found a valid sink
            }
        }
        
        // If no valid sink found with any neighbor, skip this root
        if (t == -1) {
            continue;
        }
        
        auto e1 = chrono::high_resolution_clock::now();
        ItaiRodeh ir(h, r, t, st_num);
        auto [parS, parT] = ir.construct();
        auto e2 = chrono::high_resolution_clock::now();
        eist_ms += chrono::duration<double, milli>(e2 - e1).count();
        
        double cost = compute_cost(parS, parT, r);
        
        if (cost < best_cost) {
            best_cost = cost;
            best_parS = parS;
            best_parT = parT;
            best_root = r;
            best_neighbor = t;
        }
    }
    
    if (best_root == -1) {
        cout << "No valid root found" << endl;
        return {vector<int>(), vector<int>(), -1, -1, -1.0, cand_ms, st_ms, eist_ms};
    }

    cout << "Best root: " << best_root 
         << ", neighbor: " << best_neighbor 
         << ", cost: " << best_cost << endl;
    
    return {best_parS, best_parT, best_root, best_neighbor, best_cost, cand_ms, st_ms, eist_ms};
}



int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <edge_list_file>" << endl;
        return 1;
    }
    
    string filename = argv[1];
    cout << "Opening file: " << filename << endl;
    
    ifstream infile(filename);
    
    if (!infile) {
        cerr << "Cannot open file: " << filename << endl;
        return 1;
    }
    
    // Load graph and remap IDs
    unordered_map<int, int> id_map;
    vector<pair<int, int>> edges;
    string line;
    int line_count = 0;
    
    while (getline(infile, line)) {
        line_count++;
        if (line.empty() || line[0] == '#') continue;
        int u, v;
        if (sscanf(line.c_str(), "%d %d", &u, &v) == 2) {
            if (id_map.find(u) == id_map.end()) {
                id_map[u] = id_map.size();
            }
            if (id_map.find(v) == id_map.end()) {
                id_map[v] = id_map.size();
            }
            edges.push_back({id_map[u], id_map[v]});
        }
    }
    infile.close();
    
    cout << "Read " << line_count << " lines, " << edges.size() << " edges, " 
         << id_map.size() << " vertices" << endl;
    
    int n = id_map.size();
    Graph g(n);
    for (auto& e : edges) {
        g.add_edge(e.first, e.second);
    }
    
    cout << "Building graph complete" << endl;
    
    auto total_start = chrono::high_resolution_clock::now();
    
    // Extract largest BiCC
    auto b1 = chrono::high_resolution_clock::now();
    BiCCExtractor extractor(g);
    Graph h = extractor.extract_largest();
    auto b2 = chrono::high_resolution_clock::now();
    double bicc_ms = chrono::duration<double, milli>(b2 - b1).count();
    
    cout << "Largest BiCC: " << h.n << " vertices" << endl;
    
    if (h.n == 0) {
        cerr << "No biconnected component found" << endl;
        return 1;
    }
    
    int edges_count = 0;
    for (int i = 0; i < h.n; i++) {
        edges_count += h.adj[i].size();
    }
    edges_count /= 2;
    
    cout << "BiCC edges: " << edges_count << endl;
    
// ========== ROOT-AWARE VARIANT ==========

auto result = root_aware_eist(h);

auto total_end = chrono::high_resolution_clock::now();
double total_runtime_ms = chrono::duration<double, milli>(total_end - total_start).count();

cout << "BiCC: " << bicc_ms << " ms" << endl;
cout << "Candidates: " << result.cand_ms << " ms" << endl;
cout << "ST total: " << result.st_ms << " ms" << endl;
cout << "EIST total: " << result.eist_ms << " ms" << endl;
cout << "Total root-aware runtime: " << total_runtime_ms << " ms" << endl;

// Compute depth metrics with correct root
auto [S_max, S_avg] = compute_depth_metrics(result.parS, result.root);
auto [T_max, T_avg] = compute_depth_metrics(result.parT, result.root);

cout << "S: max=" << S_max << ", avg=" << S_avg << endl;
cout << "T: max=" << T_max << ", avg=" << T_avg << endl;
    
    // Extract dataset name from filename
    string dataset = filename;
    size_t last_slash = dataset.find_last_of("/\\");
    if (last_slash != string::npos) {
        dataset = dataset.substr(last_slash + 1);
    }
    
    // Append to CSV
    ofstream csvfile("results/results_sequential.csv", ios::app);
    
    // Check if file is empty to write header
    ifstream check("results/results_sequential.csv");
    bool is_empty = check.peek() == ifstream::traits_type::eof();
    check.close();
    
    if (is_empty) {
        csvfile << "dataset,algorithm,V,E,bicc_ms,cand_ms,st_ms,eist_ms,total_runtime_ms,cost,S_max,S_avg,T_max,T_avg" << endl;
    }
    
    csvfile << dataset << ",RA," << h.n << "," << edges_count << "," 
            << fixed << setprecision(2) << bicc_ms << ","
            << result.cand_ms << "," << result.st_ms << "," << result.eist_ms << ","
            << total_runtime_ms << ","
            << setprecision(2) << result.cost << ","
            << S_max << "," << setprecision(1) << S_avg << ","
            << T_max << "," << T_avg << endl;
    
    csvfile.close();
    
    cout << "Results written to results/results_sequential.csv" << endl;
    
    return 0;
}
