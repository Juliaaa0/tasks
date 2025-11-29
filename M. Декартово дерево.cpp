#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Node {
    int n, k, pr;
    Node *l, *r, *p;

    Node() {
        k = 0;
        pr = 0;
        n = 0;
        l = nullptr;
        r = nullptr;
        p = nullptr;
    }
    Node(int z, int x, int y) {
        n = z;
        k = x;
        pr = y;
        l = nullptr;
        r = nullptr;
        p = nullptr;

    }
};

bool compare(const Node& a, const Node& b) {
    return a.k < b.k;
}

void Build(vector<Node>& V, int n, vector<vector<int>>& ans) {
    Node* prev = &V[0];
    for (int i = 1; i < n; ++i) {
            //cout << V[i].k << " " << V[i].pr << " " << V[i].n << endl;
            if (V[i].pr > prev->pr) {
                prev->r = &V[i];
                V[i].p = prev;
                ans[V[i].n][0] = prev->n;
                ans[prev->n][2] = V[i].n;
                //cout << 1 << " " << prev->n << endl;
            }
            else {
                Node *v = prev;
                while (v->p != nullptr && v->pr > V[i].pr) {
                    v = v->p;
                }
                //cout << v->n << endl;

                if (v->pr < V[i].pr) {
                    v->r->p = &V[i];
                    V[i].l = v->r;
                    v->r = &V[i];
                    V[i].p = v;
                    ans[v->n][2] = V[i].n;
                    ans[V[i].n][0] = v->n;
                    ans[V[i].l->n][0] = V[i].n;
                    ans[V[i].n][1] = V[i].l->n;
                    //cout << 2 << " " << V[i].n << endl;
                }
                else {
                    v->p = &V[i];
                    V[i].l = v;
                    ans[v->n][0] = V[i].n;
                    ans[V[i].n][1] = v->n;
                    //cout << 3 << " " << V[i].n << endl;
                }
                //cout << 00 << endl;
            }
            prev = &V[i];
    }
}

int main()
{
    int N, a, b;
    vector<Node> nodes;
    cin >> N;

    for (int i = 1; i <= N; i++) {
        cin >> a >> b;
        nodes.push_back(Node(i, a, b));
    }

    vector<vector<int>> ans(N+1, vector<int>(3, 0));
    sort(nodes.begin(), nodes.end(), compare);


    Node* T = nullptr;

    Build(nodes, N, ans);

    cout << "YES" << endl;
    for (int j = 1; j <= N; ++j) {
        cout << ans[j][0] << " " << ans[j][1] << " " << ans[j][2] << endl;
    }

    return 0;
}
