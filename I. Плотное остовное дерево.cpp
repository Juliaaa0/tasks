#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Node {
    int h = 0, r;

    Node() {
        r = 0;
    }

    Node(int u) {
        r = u;
    }
};

int Root(vector<Node>& V, int u) {
    int x;
    if (V[u].r == u) {
        return u;
    }
    else {
        x = Root(V, V[u].r);
        V[u].r = x;
        return x;
    }
}

bool Kruskal(vector<vector<int>>& ed, vector<Node>& V, int m, int n, int& d, int j) {

    for (int i = 0; i <= n; ++i){
        V[i].r = i;
        V[i].h = 0;
    }

    int b, e, w, k = 0, x, y;
    for (int i = j; i < m; ++i) {
        b = ed[i][0];
        e = ed[i][1];
        w = ed[i][2];
        x = Root(V, b);
        y = Root(V, e);
        //cout << b << " " << e << " " << x << " " << y << endl;
        if (x != y) {
            k += 1;
            if(k == n - 1) {
                d = w;
                return true;
            }

            if (V[x].h == V[y].h) {
                V[x].r = y;
                V[y].h += 1;
            }
            else if (V[x].h < V[y].h) {
                V[x].r = y;
            }
            else {
                V[y].r = x;
            }
        }
    }
    return false;
}

int main()
{
    int n, m;
    cin >> n >> m;

    vector<Node> V(n+1);
    for (int i = 1; i <= n; ++i){
        V[i] = Node(i);
    }

    int b, e, w;
    vector<vector<int>> ed(m);
    for (int i = 0; i < m; ++i) {
        cin >> b >> e >> w;
        ed[i].push_back(b);
        ed[i].push_back(e);
        ed[i].push_back(w);
    }

    sort(ed.begin(), ed.end(), [](const vector<int>& v1, const vector<int>& v2) {
        return v1[2] < v2[2];
    });

    int ans = 2000000001, t = -1000000001, d = 2000000001;

    for (int i = 0; i < m - n + 2; ++i) {

        b = ed[i][0];
        e = ed[i][1];
        w = ed[i][2];
        int f = ed[i + n - 2][2];
        //cout << w << " " << t << " " << f - w << " " << ans << endl;
        if (w == t || ((f - w) >= ans)) {
            continue;
        }
        else {
            //cout << i << endl;
            if (Kruskal(ed, V, m, n, d, i)) {
                //cout << i << " " << d << endl;
                if ((d - w) < ans) {
                    ans = d - w;
                }
            }
            t = w;
        }
    }

    if (ans < 2000000001) {
        cout << "YES" << endl << ans;
    }
    else {
        cout << "NO";
    }

    return 0;
}
