#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

bool DFS(int u, int t, vector<vector<int>>& adj, vector<int>& d, vector<int>& path) {
    d[u] = 1;
    if (u == t) {
        path.push_back(u);
        return true;
    }

    for (int i = 0; i < adj[u].size(); ++i) {
        int v = adj[u][i];
        if (d[v] == 0 && DFS(v, t, adj, d, path)) {
            path.push_back(u);
            return true;
        }
    }

    return false;
}

int main()
{
    int n, m, v, u;
    cin >> n >> m;

    vector<vector<int>> adj(n + m + 2);
    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            cin >> v;
            if (v > 0) {
                adj[i].push_back(v + n);
            }
            else {
                break;
            }
        }
    }

    for (int i = 1; i <= n; ++i) {
        adj[0].push_back(i);
    }

    int t = m + n + 1;
    for (int i = n + 1; i < t; ++i) {
        adj[i].push_back(t);
    }

    vector<int> d(t + 1);
    vector<int> path;
    int k = 1;
    while (k == 1) {
        for (int i = 0; i <= t; ++i) {
            d[i] = 0;
        }
        path.clear();

        if (DFS(0, t, adj, d, path)) {
            reverse(path.begin(), path.end());
            //for (int h = 0; h < path.size(); ++h) {
                //cout << path[h] << " ";
            //}
            //cout << endl;
            for (int j = 0; j < path.size() - 1; ++j) {
                u = path[j];
                v = path[j + 1];
                adj[u].erase(remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
                adj[v].push_back(u);
            }
        }
        else {
            k = 0;
        }
    }

    int ans = 0;
    for (int i = n + 1; i <= n + m; ++i) {
        for (int v : adj[i]) {
            if (v != t) {
                ++ans;
                break;
            }
        }
    }
    cout << ans << endl;

    for (int i = n + 1; i <= n + m; ++i) {
        for (int v : adj[i]) {
            if (v != t) {
                cout << v << " " << i - n << endl;
                break;
            }
        }
    }

    return 0;
}
