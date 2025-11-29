#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

vector<int> ver;

void DFS1(int u, vector<vector<int>>& adj, vector<int>& d) {
    //cout << u << endl;
    int v;
    d[u] = 1;
    for (int i = 0; i < adj[u].size(); ++i) {
        //cout << adj[u][i] << " " << d[adj[u][i]] << endl;
        v = adj[u][i];
        if (d[v] == 0) {
            DFS1(v, adj, d);
        }
    }
    ver.push_back(u);
}

int main()
{
    int n, m, a, b;
    cin >> n >> m;

    vector<int> d(n, 0);
    vector<vector<int>> adj(n);
    for (int i = 0; i < m; ++i) {
        cin >> a >> b;
        --a;
        --b;
        adj[a].push_back(b);
    }

    for (int i = 0; i < n; ++i) {
        if (d[i] == 0) {
            DFS1(i, adj, d);
        }
    }

    int ans = 1, u, v, k;
    for (int i = 0; i < n - 1; ++i) {
        u = ver[n - 1 - i];
        v = ver[n - 2 - i];
        k = 0;
        for (int j = 0; j < adj[u].size(); ++j) {
            if (adj[u][j] == v) {
                k = 1;
                break;
            }
        }
        if (k == 0) {
            ans = 0;
            break;
        }
    }

    if (ans == 1) {
        cout << "YES";
    }
    else {
        cout << "NO";
    }

    return 0;
}
