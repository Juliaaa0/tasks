#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <memory>
#include <bitset>

//using namespace std;

std::pair<std::vector<int>, std::vector<int>> AllocateV3(int jobspec,  int status[6][12]) {
    class Submatrix {
    public:
        std::vector<int> rows;
        std::vector<int> columns;
        Submatrix(const std::vector<int>& rows, const std::vector<int>& columns)
            : rows(rows), columns(columns) {}
        // проверка на пересечение двух подматриц
        bool intersects(const Submatrix& other) const {
            std::vector<int> common_rows, common_cols;

            set_intersection(rows.begin(), rows.end(),
                                 other.rows.begin(), other.rows.end(),
                                 back_inserter(common_rows));
            set_intersection(columns.begin(), columns.end(),
                                 other.columns.begin(), other.columns.end(),
                                 back_inserter(common_cols));

            return !common_rows.empty() && !common_cols.empty();
        }

        // содержитс€ ли в self подматрица other
        bool contains(const Submatrix& other) const {
            return includes(rows.begin(), rows.end(), other.rows.begin(), other.rows.end()) && includes(columns.begin(), columns.end(), other.columns.begin(), other.columns.end());
        }
    };

    class MaxSubmatrix {
    public:
        std::vector<Submatrix> max_rects;
        std::vector<Submatrix> jobs;
        std::vector<int> row_ins; // количество зан€тых €чеек в строках
        std::vector<int> col_ins; // количество зан€тых €чеек в столбцах
        std::vector<std::vector<int>> col_num; // номера строк, в которых заблокированы €чейки дл€ каждого столца
        std::vector<std::vector<int>> row_num; // номера столбцов, в которых заблокированы €чейки дл€ каждой строки

        MaxSubmatrix(const std::vector<Submatrix>& max_rects,
                     const std::vector<Submatrix>& jobs,
                     const std::vector<int>& row_ins,
                     const std::vector<int>& col_ins,
                     const std::vector<std::vector<int>>& col_num,
                     const std::vector<std::vector<int>>& row_num)
            : max_rects(max_rects), jobs(jobs), row_ins(row_ins), col_ins(col_ins), col_num(col_num), row_num(row_num) {}

    private:
        bool _rect_fitness(const Submatrix& max_rect, int S) {
            int s = max_rect.rows.size() * max_rect.columns.size();
            return s >= S;
        }

//        std::tuple<int, int, int> pos(const Submatrix& rect, int S) {
//            int cnt_row = rect.rows.size();
//            int cnt_col = rect.columns.size();
//
//            for (int i = 1; i <= std::min(S, cnt_row); i++) {
//                int n = ceil(static_cast<double>(S) / i);
//                if (n <= cnt_col) {
//                    return std::make_tuple(i, n, i * n - S);
//                }
//            }
//            return std::make_tuple(0, 0, 0);
//        }
         std::tuple<int, int, int> pos(const Submatrix& rect, int S) {
            int cnt_row = rect.rows.size();
            int cnt_col = rect.columns.size();
            int n0, m0, s0 = cnt_row * cnt_col;

            std::vector<std::tuple<int, int, int>> var;

            for (int i = std::min(S, cnt_row); i >= 1; --i) {
                int n = ceil(static_cast<double>(S) / i);
                if (n <= cnt_col && i * n - S < s0) {
                    n0 = i;
                    m0 = n;
                    s0 =  i * n - S;
                }
            }

            return std::make_tuple(n0, m0, s0);
        }

        Submatrix best_place(const Submatrix& rect, int n, int m) {
            std::vector<std::pair<int, int>> best_row;
            for (int row : rect.rows) {
                best_row.push_back({row, row_ins[row]});
            }

            std::vector<std::pair<int, int>> best_col;
            for (int col : rect.columns) {
                best_col.push_back({col, col_ins[col]});
            }

            sort(best_row.begin(), best_row.end(),
                     [](const auto& a, const auto& b) {
                         return a.second != b.second ? a.second > b.second : a.first < b.first;
                     });
            sort(best_col.begin(), best_col.end(),
                     [](const auto& a, const auto& b) {
                         return a.second != b.second ? a.second > b.second : a.first < b.first;
                     });

            std::vector<int> selected_rows, selected_cols;
            for (int i = 0; i < n; i++) {
                selected_rows.push_back(best_row[i].first);
            }
            for (int i = 0; i < m; i++) {
                selected_cols.push_back(best_col[i].first);
            }

            sort(selected_rows.begin(), selected_rows.end());
            sort(selected_cols.begin(), selected_cols.end());

            return Submatrix(selected_rows, selected_cols);
        }

        uint64_t vectorToMask(const std::vector<int>& vec) {
            uint64_t mask = 0;
            for (int x : vec) {
                mask |= (1ULL << x);
                std::cout << x << " " << mask << std::endl;
            }
            return mask;
        }

        std::vector<int> best_vec_col(const std::vector<std::vector<int>>& vectors, int k) {
            int n = vectors.size();
            std::cout << "vectors.size(); " << n  << std::endl;

            std::vector<uint64_t> masks(n);
            for (int i = 0; i < n; ++i) {
                masks[i] = vectorToMask(vectors[i]);
            }

            // находим начальный вектор
            int bestStart = 0;
            int maxTotalIntersection = -1;

            for (int i = 0; i < n; i++) {
                int totalIntersection = 0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
//                        std::cout << "i j " << i << " " << j <<std::endl;
//                        std::cout << "masks[i] & masks[j] " << masks[i] << " " << masks[j] <<std::endl;
//                        std::cout << " __builtin_popcountll(masks[i] & masks[j]) " <<  __builtin_popcountll(masks[i] & masks[j]) << std::endl;
                        totalIntersection += __builtin_popcountll(masks[i] & masks[j]);
                    }
                }

//                std::cout << "bestStart " <<
                if (totalIntersection > maxTotalIntersection) {
                    maxTotalIntersection = totalIntersection;
                    bestStart = i;
                }
            }

//            std::cout << "bestStart " << bestStart << std::endl;
            // добавл€ем остальные вектроры по количеству пересечений с общим
            std::vector<bool> selected(n, false);
            uint64_t currentUnion = masks[bestStart];
            std::vector<int> result = {bestStart};
            selected[bestStart] = true;

            for (int step = 1; step < k; step++) {
                int bestIdx = -1;
                int maxIntersection = -1;

                for (int i = 0; i < n; i++) {
                    if (selected[i]) continue;

                    // —читаем пересечение с текущим объединением
                    int intersection = __builtin_popcountll(masks[i] & currentUnion);

                    if (intersection > maxIntersection) {
                        maxIntersection = intersection;
                        bestIdx = i;
                    }
                }

                if (bestIdx != -1) {
                    selected[bestIdx] = true;
                    currentUnion |= masks[bestIdx];
                    result.push_back(bestIdx);
                } else {
                    // ≈сли все оставшиес€ векторы не пересекаютс€, берем любой
                    for (int i = 0; i < n; i++) {
                        if (!selected[i]) {
                            selected[i] = true;
                            currentUnion |= masks[i];
                            result.push_back(i);
                            break;
                        }
                    }
                }
            }

            return result;
        }

        std::bitset<128> vectorToBitset(const std::vector<int>& vec) {
            std::bitset<128> bs;
            for (int x : vec) {
                bs.set(x);
            }
            return bs;
        }

        std::vector<int> best_vec_row(const std::vector<std::vector<int>>& vectors, int k) {
//            std::cout << "AAAAAAAAAAAAAAAAAA" << std::endl;
            int n = vectors.size();
//            std::cout << "vectors.size(); " << n  << std::endl;

            std::vector<std::bitset<128>> masks(n);
            for (int i = 0; i < n; ++i) {
                masks[i] = vectorToBitset(vectors[i]);
            }

            // находим начальный вектор
            int bestStart = 0;
            int maxTotalIntersection = -1;

            for (int i = 0; i < n; i++) {
                int totalIntersection = 0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
//                        std::cout << "i j " << i << " " << j <<std::endl;
//                        std::cout << "masks[i] & masks[j] " << masks[i] << " " << masks[j] <<std::endl;
//                        std::cout << " __builtin_popcountll(masks[i] & masks[j]) " <<  __builtin_popcountll(masks[i] & masks[j]) << std::endl;
                        totalIntersection += (masks[i] & masks[j]).count();
                    }
                }

//                std::cout << "bestStart " <<
                if (totalIntersection > maxTotalIntersection) {
                    maxTotalIntersection = totalIntersection;
                    bestStart = i;
                }
            }

//            std::cout << "bestStart " << bestStart << std::endl;
            // добавл€ем остальные вектроры по количеству пересечений с общим
            std::vector<bool> selected(n, false);
            std::bitset<128> currentUnion = masks[bestStart];
            std::vector<int> result = {bestStart};
            selected[bestStart] = true;

            for (int step = 1; step < k; step++) {
                int bestIdx = -1;
                int maxIntersection = -1;

                for (int i = 0; i < n; i++) {
                    if (selected[i]) continue;

                    // —читаем пересечение с текущим объединением
                    int intersection = (masks[i] & currentUnion).count();

                    if (intersection > maxIntersection) {
                        maxIntersection = intersection;
                        bestIdx = i;
                    }
                }

                if (bestIdx != -1) {
                    selected[bestIdx] = true;
                    currentUnion |= masks[bestIdx];
                    result.push_back(bestIdx);
                } else {
                    // ≈сли все оставшиес€ векторы не пересекаютс€, берем любой
                    for (int i = 0; i < n; i++) {
                        if (!selected[i]) {
                            selected[i] = true;
                            currentUnion |= masks[i];
                            result.push_back(i);
                            break;
                        }
                    }
                }
            }

            return result;
        }

         // теперь смотрю не на количество заблокированных €чеек в столбцах, а на количество общих строк у этих заблокированных €чеек
        Submatrix best_place2(const Submatrix& rect, int n, int m) {
            std::cout << "Submatrix& rect " << rect.rows.size() << " " << rect.columns.size() << std::endl;
            std::cout << "sizes " << n << " " << m << std::endl;
            std::vector<int> selected_rows, selected_cols;

            std::vector<std::pair<int, int>> all_col;
            for (int col : rect.columns) {
                all_col.push_back({col, col_num[col].size()});
            }
            std::vector<std::pair<int, int>> all_row;
            for (int row : rect.rows) {
                all_row.push_back({row, row_num[row].size()});
            }

            auto col2vec = make2vectors(all_col, m);
//            std::cout << "m " << m << std::endl;
            std::cout << "col2vec.sizes " << col2vec.first.size() << " " << col2vec.second.size() << std::endl;
            auto row2vec = make2vectors(all_row, n);
            std::cout << "row2vec.sizes " << row2vec.first.size() << " " << row2vec.second.size() << std::endl;

//            std::cout << "AAA" << std::endl;
            if (col2vec.first.size() + col2vec.second.size() == m) {
                for (auto c : col2vec.first) {
                    selected_cols.push_back(c.first);
                }
                for (auto c : col2vec.second) {
                    selected_cols.push_back(c.first);
                }
            }
            else {

                for (auto c : col2vec.first) {
                    selected_cols.push_back(c.first);
                }

                int k = m - col2vec.first.size();

//                std::cout << "col2vec.second" << std::endl;
//                for (auto i : col2vec.second) {
//                    std::cout << i.first << " ";
//                }
//                std::cout << std::endl;

                std::vector<std::vector<int>> equal_cols;
                for (auto c : col2vec.second) {
                    equal_cols.push_back(col_num[c.first]);
//                    std::cout << c.first << std::endl;
                }

//                std::cout << "AAAA " << col2vec.second.size() << std::endl;
                auto col_n = best_vec_col(equal_cols, k);
//                std::cout << "AAAAA" << std::endl;
                std::sort(col_n.begin(), col_n.end());
                for (auto& c : col_n) {
                    selected_cols.push_back(col2vec.second[c].first);
//                    std::cout << col2vec.second[c].first << std::endl;
                }
            }

//
//            std::cout << "AAAAA" <<std::endl;
//            std::cout << "row2vec " << std::endl;
            if (row2vec.first.size() + row2vec.second.size() == n) {
//                std::cout << "AAAAAA" <<std::endl;
                for (auto r : row2vec.first) {
                    selected_rows.push_back(r.first);
                }
                for (auto r : row2vec.second) {
                    selected_rows.push_back(r.first);
                }
            }
            else {
                std::cout << "AAAAAA" << std::endl;
                for (auto r : row2vec.first) {
                    selected_rows.push_back(r.first);
                }

                int k = n - row2vec.first.size();

                std::vector<std::vector<int>> equal_rows;
                for (auto r : row2vec.second) {
                    equal_rows.push_back(row_num[r.first]);
                }

                auto row_n = best_vec_row(equal_rows, k);
                std::sort(row_n.begin(), row_n.end());
                for (auto& r : row_n) {
                    selected_rows.push_back(row2vec.second[r].first);
                }
            }

            return Submatrix(selected_rows, selected_cols);
        }

        std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> make2vectors(std::vector<std::pair<int, int>> pairs, int k) {
            std::vector<std::pair<int, int>> greater_pairs, equal_pairs;
            std::sort(pairs.begin(), pairs.end(),
              [](const auto& a, const auto& b) {
                  return a.second < b.second;
              });

            for (auto& i : pairs) {
                std::cout << "i " << i.first << " i_ins " << i.second << std::endl;
            }

            int k_max = pairs[k-1].second;
            std::cout << "k_max " << k_max << " k " << k << std::endl;
            for (auto& p : pairs) {
                if (p.second < k_max) {
                    greater_pairs.push_back(p);
                }
                else if (p.second == k_max) {
                    equal_pairs.push_back(p);
                }
                else {
                    break;
                }
            }

            return make_pair(greater_pairs, equal_pairs);
        }

        // теперь смотрю не на количество заблокированных €чеек в столбцах, а на количество общих строк у этих заблокированных €чеек
//        Submatrix best_place2(const Submatrix& rect, int n, int m) {
//            std::vector<std::pair<int, int>> best_row;
//            for (int row : rect.rows) {
//                best_row.push_back({row, row_ins[row]});
//            }
//
//            std::vector<std::pair<int, int>> best_col;
//            for (int col : rect.columns) {
//                best_col.push_back({col, col_ins[col]});
//            }
//
//            sort(best_row.begin(), best_row.end(),
//                     [](const auto& a, const auto& b) {
//                         return a.second != b.second ? a.second > b.second : a.first < b.first;
//                     });
//            sort(best_col.begin(), best_col.end(),
//                     [](const auto& a, const auto& b) {
//                         return a.second != b.second ? a.second > b.second : a.first < b.first;
//                     });
//
//            int cnt_cells = best_col[n-1].second;
//            std::vector<int> selected_rows, selected_cols;
//            for (int i = 0; i < n; i++) {
//                selected_rows.push_back(best_row[i].first);
//            }
//
//            std::vector<std::pair<int, int>> result;
//            auto lower = std::lower_bound(arr.begin(), arr.end(), m,
//            [](const std::pair<int, int>& p, int value) {
//                return p.second > value;
//            });
//
//            auto upper = std::upper_bound(arr.begin(), arr.end(), m,
//            [](int value, const std::pair<int, int>& p) {
//                return value > p.second;
//            });
//            result.insert(result.end(), lower, upper);
//
//            sort(selected_rows.begin(), selected_rows.end());
//            sort(selected_cols.begin(), selected_cols.end());
//
//            return Submatrix(selected_rows, selected_cols);
//        }

    public:
        std::vector<Submatrix> _remove_duplicates(std::vector<Submatrix> rects) {

            std::vector<Submatrix> new_rects;
            std::vector<bool> contained(max_rects.size(), false);
            for (size_t i = 0; i < rects.size(); ++i) {
                if (contained[i]) continue;
                for (size_t j = i + 1; j < rects.size(); ++j) {
                    if (contained[j]) continue;
                    if (rects[i].contains(rects[j])) {
                        contained[j] = true;
                    }
                    else if (rects[j].contains(rects[i])) {
                        contained[i] = true;
                    }
                }
            }

            for (int i = 0; i <  rects.size(); ++i) {
                if (!contained[i]) {
                    new_rects.push_back(rects[i]);
                }
            }

            return new_rects;
        }

//        Submatrix _select_position4(int S) {
//            if (max_rects.empty()) {
//                return Submatrix({}, {});
//            }
//            std::vector<std::pair<Submatrix, double>> fit;
//            for (const auto& m : max_rects) {
//                int s = m.rows.size() * m.columns.size();
//                if (s >= S) {
//                    double k = static_cast<double>(std::max(m.columns.size(), m.rows.size())) / std::min(m.columns.size(), m.rows.size());
//                    fit.push_back({m, k});
//                }
//            }
//
//            if (fit.empty()) {
//                return Submatrix({}, {});
//            }
//
//            sort(fit.begin(), fit.end(),
//                     [](const auto& a, const auto& b) { return a.second < b.second; });
//
//            const Submatrix& rec1 = fit[0].first;
//            const Submatrix& rec2 = fit.back().first;
//
//            auto [n1, m1, c1] = pos(rec1, S);
//            auto [n2, m2, c2] = pos(rec2, S);
//
//            if (c1 < c2) {
//                return best_place(rec1, n1, m1);
//            }
//            return best_place(rec2, n2, m2);
//

        Submatrix _select_position4(int S) {
            if (max_rects.empty()) {
                return Submatrix({}, {});
            }

            std::vector<std::pair<Submatrix, double>> fit;

            for (const auto& m : max_rects) {
                int s = m.rows.size() * m.columns.size();
                if (s >= S) {
                    double k = static_cast<double>(std::max(m.columns.size(), m.rows.size())) / std::min(m.columns.size(), m.rows.size());
                    fit.push_back({m, k});
                }
            }

            if (fit.empty()) {
                return Submatrix({}, {});
            }

            sort(fit.begin(), fit.end(),
                        [](const auto& a, const auto& b) { return a.second < b.second; });

            const Submatrix& rec1 = fit[0].first;
            const Submatrix& rec2 = fit.back().first;

            auto [n1, m1, c1] = pos(rec1, S);
            auto [n2, m2, c2] = pos(rec2, S);

            if (c1 < c2) {
                return best_place2(rec1, n1, m1);
            }
            return best_place2(rec2, n2, m2);
        }

        Submatrix _select_position_simple(int S) {
            if (max_rects.empty()) {
                return Submatrix({}, {});
            }
            for (const auto& mr : max_rects) {
                int s = mr.rows.size() * mr.columns.size();
                if (s >= S) {
                    auto [n, m, c] = pos(mr, S);
                    std::vector<int> row(mr.rows.begin(), mr.rows.begin() + n);
                    std::vector<int> col(mr.columns.begin(), mr.columns.begin() + m);
                    return Submatrix(row, col);
                }
            }
            return Submatrix({}, {});
        }

        std::pair<Submatrix, Submatrix> _generate_splits(const Submatrix& m, const Submatrix& r) {
            Submatrix r1({}, {});
            Submatrix r2({}, {});
            std::vector<int> diff_rows, diff_cols;
            set_difference(m.rows.begin(), m.rows.end(), r.rows.begin(), r.rows.end(), back_inserter(diff_rows));
            set_difference(m.columns.begin(), m.columns.end(), r.columns.begin(), r.columns.end(), back_inserter(diff_cols));
            if (!diff_rows.empty()) {
                r1 = Submatrix(diff_rows, m.columns);
            }

            if (!diff_cols.empty()) {
                r2 = Submatrix(m.rows, diff_cols);
            }

            auto p = std::make_pair(r1, r2);
            return p;
        }

//        void _split(const Submatrix& rect) {
//            std::vector<Submatrix> new_max_rects;
//
//            for (const auto& r : max_rects) {
//                if (r.intersects(rect)) {
//                    auto splits = _generate_splits(r, rect);
//                    new_max_rects.insert(new_max_rects.end(), splits.begin(), splits.end());
//                }
//                else {
//                    new_max_rects.push_back(r);
//                }
//            }
//
//            max_rects = new_max_rects;
//        }

        void _split(const Submatrix& rect) {
            std::vector<Submatrix> new_max_rects;
            std::vector<Submatrix> old_max_rects;
            std::vector<Submatrix> new_r1;
            std::vector<Submatrix> new_r2;
            int k;

            for (const auto& r : max_rects) {
                if (r.intersects(rect)) {
                    auto rects = _generate_splits(r, rect);
                    auto r1 = rects.first;
                    auto r2 = rects.second;
                    if (!r1.rows.empty()) {
                        new_r1.push_back(r1);
                    }
                    if (!r2.rows.empty()) {
                        new_r2.push_back(r2);
                    }
//                    new_max_rects.insert(new_max_rects.end(), splits.begin(), splits.end());
                }
                else {
                    old_max_rects.push_back(r);
                    new_max_rects.push_back(r);
                }
            }

            auto new_r1_1 = _remove_duplicates(new_r1);
            auto new_r2_1 = _remove_duplicates(new_r2);
            for (Submatrix r : new_r1_1) {
                k = 1;
                for (Submatrix m : old_max_rects) {
                    if (m.contains(r)){
                         k = 0;
                         break;
                    }
                }
                if (k == 1) {
                    new_max_rects.push_back(r);
                }
            }
            for (Submatrix r : new_r2_1) {
                k = 1;
                for (Submatrix m : old_max_rects) {
                    if (m.contains(r)){
                         k = 0;
                         break;
                    }
                }
                if (k == 1) {
                    new_max_rects.push_back(r);
                }
            }

            max_rects = new_max_rects;
        }

        bool add_rect(int S, const std::string& n) {
            Submatrix rect({}, {});

            if (n == "MAX_COMPACT") {
                rect = _select_position_simple(S);
            }

            if (rect.rows.empty() || rect.columns.empty()) {
                return false;
            }

            for (int r : rect.rows) {
                row_ins[r] += rect.columns.size();
            }
            for (int c : rect.columns) {
                col_ins[c] += rect.rows.size();
            }

            jobs.push_back(rect);
            _split(rect);
//            _remove_duplicates();
            return true;
        }
    };

    const int M = 6;
    const int N = 12;

    std::vector<int> all_rows(M);
    std::vector<int> all_cols(N);
    for (int i = 0; i < M; i++) all_rows[i] = i;
    for (int i = 0; i < N; i++) all_cols[i] = i;

    std::vector<int> row_ins(M, 0);
    std::vector<int> col_ins(N, 0);
    std::vector<std::vector<int>> col_num(N);
    std::vector<std::vector<int>> row_num(M);


    Submatrix x(all_rows, all_cols);
    MaxSubmatrix MS({x}, {}, row_ins, col_ins, col_num, row_num);

    // вставл€ю все заблокированные €чейки по строкам
    for (int i : all_rows) {
        std::vector<int> bl_cells, free_cells;
        for (int j : all_cols) {
            if (status[i][j] == 1 || status[i][j] == 2) {
                bl_cells.push_back(j);
                MS.col_ins[j] += 1;
            }
            else {
                MS.col_num[j].push_back(i);
                free_cells.push_back(j);
            }
        }
        std::vector<int> row_i = {i};
        MS.row_ins[i] += bl_cells.size();
        MS.row_num[i] = free_cells;
        Submatrix y(row_i, bl_cells);
        MS._split(y);
//        MS._remove_duplicates();
    }

    //вывод всех макс подматриц
//    std::cout << "MS.max_rects " << MS.max_rects.size() << std::endl;
//    for (auto& ms : MS.max_rects) {
//        for (int& r : ms.rows) {
//            std::cout << r << " ";
//        }
//        std::cout << std::endl;
//        for (int& r : ms.columns) {
//            std::cout << r << " ";
//        }
//        std::cout << std::endl;
//        std::cout << std::endl;
//    }
//    std::cout << "MS.col_num " << MS.col_num.size() << std::endl;
//    for (auto& col : MS.col_num) {
//        for (auto& i : col) {
//            std::cout << i << " ";
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << "MS.row_num " << MS.row_num.size() << std::endl;
//    for (auto& col : MS.row_num) {
//        for (auto& i : col) {
//            std::cout << i << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "MS.row_ins " << MS.row_ins.size() << std::endl;
//    for (int i : MS.row_ins) {
//        std::cout << i << " ";
//    }
//    std::cout << std::endl;
//    std::cout << "MS.col_ins " << MS.col_ins.size() << std::endl;
//    for (int i : MS.col_ins) {
//        std::cout << i << " ";
//    }
//    std::cout << std::endl;


    std::vector<int> ans_rows;
    std::vector<int> ans_cols;
    if (MS.add_rect(jobspec, "MAX_COMPACT")) {
        ans_rows = MS.jobs[0].rows;
        ans_cols = MS.jobs[0].columns;
        return {ans_rows, ans_cols};
    }
    return {ans_rows, ans_cols};
}


int main() {
    int status[6][12] = {
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        {3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1},
        {3, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1},
        {3, 1, 3, 1, 1, 3, 3, 1, 1, 1, 1, 1},
        {3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
    };
//    int status[6][12] = {
//        {1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//        {3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1},
//        {3, 1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1},
//        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
//        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
//    };
//    auto ans = AllocateV3(10, status);
//    std::cout << ans.first.size() << " " << ans.second.size() << std::endl;

    std::vector<std::vector<char>> cl(6, std::vector<char>(12, 'O'));

    int k = 1;
    std::vector<int> jobs = {8};
    for (int job : jobs) {
        auto ans = AllocateV3(job, status);
        std::vector<int> row = ans.first;
        std::vector<int> col = ans.second;
        for (int i : row) {
            for (int j : col) {
                cl[i][j] = '0' + k;
                status[i][j] = 1;
//                std::cout << i << " " << j << std::endl;
            }
        }
        ++k;

    }
    for (const auto& row : cl) {
        for (char cell : row) {
           std:: cout << cell << " ";
        }
        std::cout << std::endl;
    }
}
