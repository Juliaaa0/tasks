#include <iostream>
#include <vector>
#include <map>
#include <set>

using namespace std;

bool exist(vector<int>& A, vector<int>& B, int d) {
    int n = A.size();
    multiset<int> ms1, ms2;
    map<multiset<int>, int> cnt;
    for (int i = 0; i < d; ++i) {
        ms1.insert(A[i]);
    }
    cnt[ms1] = 1;

    for (int i = d; i < n; ++i) {
        ms1.insert(A[i]);
        ms1.erase(ms1.find(A[i - d]));
        cnt[ms1] = 1;
    }

    int m = B.size();
    for (int i = 0; i < d; ++i) {
        ms2.insert(B[i]);
    }
    if (cnt[ms2] > 0) {
        return true;
    }

    for (int i = d; i < m; ++i) {
        ms2.insert(B[i]);
        ms2.erase(ms2.find(B[i - d]));
        if (cnt[ms2] > 0) {
            return true;
        }
    }
    return false;
}
int main()
{
    int n, m;
    cin >> n;
    vector<int> A(n);
    for (int i = 0; i < n; ++i) {
        cin >> A[i];
    }

    cin >> m;
    vector<int> B(m);
    for (int i = 0; i < m; ++i) {
        cin >> B[i];
    }

    int l = min(n, m);
    int left = 1, right = l, ans = 0;
    while (left <= right) {
        int mid = (left + right)/2;
        if (exist(A, B, mid)) {
            left = mid + 1;
            ans = mid;
        }
        else {
            right = mid - 1;
        }
    }

    cout << ans;
    return 0;
}
