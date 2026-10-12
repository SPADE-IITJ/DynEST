/*
 * Test for st-numbering using OGDF-style pathfinder algorithm
 * (Even-Tarjan as implemented in the Open Graph Drawing Framework)
 */
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include <stack>
#include <set>
#include <climits>

using namespace std;

struct Graph {
    int N;
    vector<vector<int>> adj;
    Graph(int n) : N(n), adj(n) {}
    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
};

// ---- OGDF-style st-numbering (pathfinder) ----
//
// Based on the Even-Tarjan algorithm as implemented in OGDF's STNumbering.cpp.
//
// 1. Pre-number t as dfn=1. Then DFS from s (giving s dfn=2, etc.).
//    During DFS, compute low[v], dfsParent[v], and followLowPath[v].
// 2. Mark s and t, push t then s onto a stack.
// 3. Pop v. Try to find an st-path from v by scanning v's edges:
//    - Tree edge to child w: follow lowpath from w to a marked vertex.
//    - Back edge to descendant w (dfn[v] < dfn[w]): follow tree-parent
//      edges from w up to a marked vertex.
//    If path found, push path vertices (so v is on top of stack).
//    If no path found, number v = count++.
// 4. After loop, number t = count.
vector<int> st_number(const Graph& G, int n, int s, int t) {
    if (n <= 0) return {};
    if (n == 1) return {1};
    if (n == 2) { vector<int> g(2,0); g[s]=1; g[t]=2; return g; }

    vector<int> dfn(n, 0);   // 0 = unvisited; 1-indexed DFS numbers
    vector<int> low(n, 0);
    vector<int> dfsParent(n, -1);     // parent in DFS tree
    vector<int> followLow(n, -1);     // next vertex on path toward low ancestor

    int cnt = 1;

    // Step 1a: pre-number t
    dfn[t] = cnt++;  // dfn[t] = 1
    low[t] = dfn[t];

    // Step 1b: iterative DFS from s 
    // (t is already visited, so edge s-t becomes a back edge)
    struct Frame { int u; int idx; };
    stack<Frame> dfs_stk;
    dfn[s] = cnt++;  // dfn[s] = 2
    low[s] = dfn[s];
    dfs_stk.push({s, 0});

    while (!dfs_stk.empty()) {
        auto& fr = dfs_stk.top();
        int u = fr.u;
        bool pushed = false;

        while (fr.idx < (int)G.adj[u].size()) {
            int w = G.adj[u][fr.idx];
            fr.idx++;
            if (dfn[w] == 0) {
                // Tree edge
                dfsParent[w] = u;
                dfn[w] = cnt++;
                low[w] = dfn[w];
                dfs_stk.push({w, 0});
                pushed = true;
                break;
            } else if (low[u] > dfn[w]) {
                low[u] = dfn[w];
                followLow[u] = w;
            }
        }

        if (!pushed) {
            dfs_stk.pop();
            // Update parent's low
            int p = dfsParent[u];
            if (p >= 0 && low[p] > low[u]) {
                low[p] = low[u];
                followLow[p] = u;
            }
        }
    }

    // Ensure low[t] is correct
    if (low[t] > low[s]) low[t] = low[s];

    // Step 2: Pathfinder phase
    vector<bool> markedNode(n, false);
    markedNode[s] = markedNode[t] = true;

    // Edge marking: use set of canonical edge keys
    set<pair<int,int>> markedEdgeSet;
    auto markE = [&](int a, int b) {
        markedEdgeSet.insert({min(a,b), max(a,b)});
    };
    auto isMarkedE = [&](int a, int b) -> bool {
        return markedEdgeSet.count({min(a,b), max(a,b)});
    };

    markE(s, t);  // mark the st edge

    // Per-vertex adjacency iterator (to resume scanning)
    vector<int> adjIdx(n, 0);

    stack<int> nodeStack;
    nodeStack.push(t);  // bottom
    nodeStack.push(s);  // top

    vector<int> g(n, 0);
    int num = 1;

    int v = nodeStack.top(); nodeStack.pop();

    while (v != t) {
        bool pathFound = false;

        while (adjIdx[v] < (int)G.adj[v].size()) {
            int w = G.adj[v][adjIdx[v]];
            adjIdx[v]++;

            if (isMarkedE(v, w)) continue;
            markE(v, w);

            if (dfsParent[w] == v) {
                // Tree edge v -> child w
                // Follow lowpath from w down to a marked vertex
                vector<int> path;
                path.push_back(v);

                while (!markedNode[w]) {
                    path.push_back(w);
                    markedNode[w] = true;

                    int nxt = followLow[w];
                    if (nxt < 0) break;
                    markE(w, nxt);
                    w = nxt;
                }

                // Push path onto nodeStack in reverse so v is on top
                for (int i = (int)path.size()-1; i >= 0; i--)
                    nodeStack.push(path[i]);

                pathFound = true;
                break;

            } else if (dfn[v] < dfn[w]) {
                // Back edge: v (ancestor, lower dfn) to w (descendant, higher dfn)
                // Follow tree-parent edges from w up to a marked vertex
                vector<int> path;
                path.push_back(v);

                while (!markedNode[w]) {
                    path.push_back(w);
                    markedNode[w] = true;

                    int par = dfsParent[w];
                    if (par < 0) break;
                    markE(w, par);
                    w = par;
                }

                for (int i = (int)path.size()-1; i >= 0; i--)
                    nodeStack.push(path[i]);

                pathFound = true;
                break;
            }
            // else: back edge to ancestor (dfn[v] > dfn[w]) - skip
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

// ---- Tests ----

bool check_st(const Graph& G, const vector<int>& g, int n, int s, int t, const string& name) {
    cout << "\n" << name << " st-numbering (s=" << s << ", t=" << t << "):" << endl;
    for (int i = 0; i < n; i++) cout << "  g[" << i << "] = " << g[i] << endl;

    bool ok = true;

    if (g[s] != 1) { cout << "  FAIL: g[s]=" << g[s] << " != 1" << endl; ok = false; }
    if (g[t] != n) { cout << "  FAIL: g[t]=" << g[t] << " != " << n << endl; ok = false; }

    vector<bool> used(n+1, false);
    for (int i = 0; i < n; i++) {
        if (g[i] < 1 || g[i] > n) { cout << "  FAIL: g[" << i << "]=" << g[i] << " out of range" << endl; ok = false; }
        else if (used[g[i]]) { cout << "  FAIL: g[" << i << "]=" << g[i] << " duplicate" << endl; ok = false; }
        else used[g[i]] = true;
    }

    int failures = 0;
    for (int v = 0; v < n; v++) {
        bool hasLower = false, hasHigher = false;
        for (int w : G.adj[v]) {
            if (g[w] < g[v]) hasLower = true;
            if (g[w] > g[v]) hasHigher = true;
        }
        if (g[v] != 1 && !hasLower) { cout << "  FAIL: v=" << v << " (g=" << g[v] << ") no lower neighbor" << endl; ok = false; failures++; }
        if (g[v] != n && !hasHigher) { cout << "  FAIL: v=" << v << " (g=" << g[v] << ") no higher neighbor" << endl; ok = false; failures++; }
    }

    if (ok) cout << "  ALL CHECKS PASS!" << endl;
    else cout << "  " << failures << " st-property failures" << endl;
    return ok;
}

void test_K4() {
    Graph G(4);
    G.add_edge(0,1); G.add_edge(0,2); G.add_edge(0,3);
    G.add_edge(1,2); G.add_edge(1,3); G.add_edge(2,3);
    auto g = st_number(G, 4, 0, 1);
    check_st(G, g, 4, 0, 1, "K4");
}

void test_C5() {
    Graph G(5);
    G.add_edge(0,1); G.add_edge(1,2); G.add_edge(2,3);
    G.add_edge(3,4); G.add_edge(4,0);
    auto g = st_number(G, 5, 0, 1);
    check_st(G, g, 5, 0, 1, "C5");
}

void test_petersen() {
    Graph G(10);
    G.add_edge(0,1); G.add_edge(1,2); G.add_edge(2,3); G.add_edge(3,4); G.add_edge(4,0);
    G.add_edge(5,7); G.add_edge(7,9); G.add_edge(9,6); G.add_edge(6,8); G.add_edge(8,5);
    G.add_edge(0,5); G.add_edge(1,6); G.add_edge(2,7); G.add_edge(3,8); G.add_edge(4,9);
    auto g = st_number(G, 10, 0, 1);
    check_st(G, g, 10, 0, 1, "Petersen");
}

void test_K33() {
    Graph G(6);
    for (int i = 0; i < 3; i++)
        for (int j = 3; j < 6; j++)
            G.add_edge(i, j);
    auto g = st_number(G, 6, 0, 3);
    check_st(G, g, 6, 0, 3, "K3,3");
}

void test_wheel() {
    Graph G(7);
    for (int i = 1; i <= 6; i++) {
        G.add_edge(0, i);
        G.add_edge(i, (i%6)+1);
    }
    auto g = st_number(G, 7, 0, 1);
    check_st(G, g, 7, 0, 1, "Wheel W7");
}

int main() {
    test_K4();
    test_C5();
    test_petersen();
    test_K33();
    test_wheel();
    cout << "\n=== All tests complete ===" << endl;
    return 0;
}
