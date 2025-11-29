#include <iostream>
#include <vector>
#include <string>
#include <set>
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
        return x;
    }
}

void Merge(vector<Node>& V, int u, int v) {
    int x, y;
    x = Root(V, u);
    y = Root(V, v);
    if (x != y) {
        if (V[x].h == V[y].h) {
            V[y].r = x;
            V[y].h += 1;
        }
        else if (V[x].h > V[y].h) {
            V[y].r = x;
        }
        else {
            V[x].r = y;
        }
    }
}

bool Ans(vector<Node>& V, int x, int y) {
    if (Root(V, x) == Root(V, y)) {
        return true;
    }
    return false;
}

int main()
{
    int n, q, type, x, y;
    cin >> n >> q;

    vector<Node> V(n + 1);
    for (int i = 0; i <= n; ++i){
        V[i] = Node(i);
    }

    set<pair<int, int>> Set;
    int a, b, c;
    vector<string> ans;
    vector<pair<int, int>> com;

    for (int i = 0; i < q; ++i) {
        cin >> type >> x >> y;
        if (type == 1) {
            Merge(V, x, y);
        }
        else if (type == 2) {
            for (const auto& p : Set) {
                a = p.first;
                b = p.second;
                if ((x >= a && x <= b) || (y >= a) && (y <= b)) {
                    c = com[p];
                    p = {min(a, x), max{b, y}};
                    com[p] = c;
                }
            }
        }
        else {
            if (Ans(V, x, y)) {
                ans.push_back("YES");
            }
            else {
                ans.push_back("NO");
            }
        }
    }

    for (int i = 0; i < ans.size(); ++i) {
        cout << ans[i] << endl;
    }

    return 0;
}
