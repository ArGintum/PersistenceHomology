/*
 Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes
 MIT License
 Copyright (c) 2015–2021 Ulrich Bauer
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 You are under no obligation whatsoever to provide any bug fixes, patches, or
 upgrades to the features, functionality or performance of the source code
 ("Enhancements") to anyone; however, if you choose to make your Enhancements
 available either publicly, or directly to the author of this software, without
 imposing a separate written license agreement for such Enhancements, then you
 hereby grant the following license: a non-exclusive, royalty-free perpetual
 license to install, use, modify, prepare derivative works, incorporate into
 other computer software, distribute, and sublicense such enhancements or
 derivative works thereof, in binary and source code form.
*/



#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <thread>
#include <unordered_map>

#include <tbb/parallel_for.h>
#include <boost/dynamic_bitset.hpp>

#include <tbb/concurrent_unordered_map.h>
#include <tbb/atomic.h>
#include <tbb/parallel_sort.h>


using value_t = long double;
typedef int64_t index_t;

using entry_hash_map = tbb::concurrent_unordered_map<index_t, std::pair<tbb::atomic<index_t>, value_t>>;

typedef std::vector<boost::dynamic_bitset<> > BitMatrix;

struct diameter_index_t : std::pair<value_t, index_t> {
	using std::pair<value_t, index_t>::pair;
};

using SparseColumn = std::vector<diameter_index_t>;
using PColumn = std::shared_ptr<SparseColumn>;
using AsyncMatrix = std::vector<PColumn>;

value_t get_diameter(const diameter_index_t& i) {
    return i.first;
}

index_t get_index(const diameter_index_t& i) {
    return i.second;
}

struct greater_diameter_or_smaller_index {
	bool operator()(const diameter_index_t& a, const diameter_index_t& b) {
		return (get_diameter(a) > get_diameter(b)) ||
		       ((get_diameter(a) == get_diameter(b)) && (get_index(a) < get_index(b)));
	}
};

bool compare_diameter_entries(const diameter_index_t& a, const diameter_index_t& b) {
	return (get_diameter(a) > get_diameter(b)) ||
	       ((get_diameter(a) == get_diameter(b)) && (get_index(a) < get_index(b)));
}

void check_overflow(index_t i) {
	if (i < 0)
		throw std::overflow_error("Overflow occured");
}

class compressed_lower_distance_matrix {
public:
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	void init_rows() {
		value_t* pointer = &distances[0];
		for (index_t i = 1; i < (int)size(); ++i) {
			rows[i] = pointer;
			pointer += i;
		}
	}

    compressed_lower_distance_matrix() {};

	compressed_lower_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	value_t operator()(const index_t i, const index_t j) const {
		return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
	}

	size_t size() const { return rows.size(); }
};

class BinominalHash {
	std::vector<std::vector<index_t>> B;

public:
	BinominalHash(index_t n, index_t k) {
        B.resize(n + 1);
        for (index_t i = 0; i <= n; ++i) {
            B[i].resize(k + 1, 0);
            B[i][0] = 1;
            for (index_t j = 1; j < std::min(i, k + 1); ++j)
                B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
            if (i <= k)
                B[i][i] = 1;
        check_overflow(B[i][std::min(i >> 1, k)]);
        }
    }

	index_t operator()(index_t n, index_t k) const {
        assert(n < (int)B.size() && k < (int)B[n].size());
        return B[n][k];
    }
};


class union_find {
	std::vector<index_t> parent;
	std::vector<uint8_t> rank;

public:
	union_find(const index_t n) : parent(n), rank(n, 0) {
		for (index_t i = 0; i < n; ++i) parent[i] = i;
	}

	index_t find(index_t x) {
		index_t y = x, z;
		while ((z = parent[y]) != y)
            y = z;
		while ((z = parent[x]) != y) {
			parent[x] = y;
			x = z;
		}
		return z;
	}

	void link(index_t x, index_t y) {
		if ((x = find(x)) == (y = find(y)))
            return;
		if (rank[x] > rank[y])
			parent[y] = x;
		else {
			parent[x] = y;
			if (rank[x] == rank[y])
                ++rank[y];
		}
	}
};

template <typename Heap> diameter_index_t pop_pivot(Heap& column) {
	if (column.empty()) return diameter_index_t(0, -1);

	auto pivot = column.top();
	column.pop();
	while (!column.empty() && get_index(column.top()) == get_index(pivot)) {
		column.pop();
		if (column.empty())
			return diameter_index_t(0, -1);
		else {
			pivot = column.top();
			column.pop();
		}
	}
	return pivot;
}

template <typename Heap> diameter_index_t get_pivot(Heap& column) {
	diameter_index_t result = pop_pivot(column);
	if (get_index(result) != -1)
        column.push(result);
	return result;
}

diameter_index_t low(const SparseColumn* c)
{
    if (!c)
        return diameter_index_t(0, -1);
    return c->empty() ? diameter_index_t(0, -1) : c->back();
}

bool is_zero(const SparseColumn* c)
{
    if (!c)
        return false;

    return c->empty();
}

template <class Predicate>
index_t get_max(index_t top, const index_t bottom, const Predicate pred) {
	if (!pred(top)) {
		index_t count = top - bottom;
		while (count > 0) {
			index_t step = count >> 1, mid = top - step;
			if (!pred(mid)) {
				top = mid - 1;
				count -= step + 1;
			} else
				count = step;
		}
	}
	return top;
}

class processor {
	compressed_lower_distance_matrix dist;
	index_t n, dim_max;
	value_t threshold;
	const BinominalHash binomial_coeff;

public:
	processor(compressed_lower_distance_matrix&& _dist, index_t _dim_max, value_t _threshold)
	    : dist(std::move(_dist)), n(dist.size()),
	      dim_max(std::min(_dim_max, index_t(dist.size() - 2))), threshold(_threshold),
	      binomial_coeff(n, dim_max + 2) {}

	index_t get_max_vertex(const index_t idx, const index_t k, const index_t n) const {
		return get_max(n, k - 1, [&](index_t w) -> bool { return (binomial_coeff(w, k) <= idx); });
	}

	index_t get_edge_index(const index_t i, const index_t j) const {
		return binomial_coeff(i, 2) + j;
	}

	template <typename OutputIterator>
	OutputIterator get_simplex_vertices(index_t idx, const index_t dim, index_t n,
	                                    OutputIterator out) const {
		--n;
		for (index_t k = dim + 1; k > 0; --k) {
			n = get_max_vertex(idx, k, n);
			*out++ = n;
			idx -= binomial_coeff(n, k);
		}
		return out;
	}

	class simplex_coboundary_enumerator {
		index_t idx_below, idx_above, v, k;
		std::vector<index_t> vertices;
		const diameter_index_t simplex;
		const compressed_lower_distance_matrix& dist;
		const BinominalHash& binomial_coeff;

	public:
		simplex_coboundary_enumerator(const diameter_index_t _simplex, const index_t _dim,
		                              const processor& parent)
		    : idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1), k(_dim + 1),
		      vertices(_dim + 1), simplex(_simplex), dist(parent.dist),
		      binomial_coeff(parent.binomial_coeff) {
			parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n, vertices.rbegin());
		}

		bool has_next(bool all_cofacets = true) {
			while ((v != -1) && (binomial_coeff(v, k) <= idx_below)) {
				if (!all_cofacets) return false;
				idx_below -= binomial_coeff(v, k);
				idx_above += binomial_coeff(v, k + 1);
				--v;
				--k;
				assert(k != -1);
			}
			return v != -1;
		}

		diameter_index_t next() {
			value_t cofacet_diameter = get_diameter(simplex);
			for (index_t w : vertices) cofacet_diameter = std::max(cofacet_diameter, dist(v, w));
			index_t cofacet_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
			return diameter_index_t(cofacet_diameter, cofacet_index);
		}
	};

	void parallel_assembling(std::vector<diameter_index_t>& simplices, std::atomic<int>& next_simplex, std::vector<diameter_index_t>& outp,
                                std::vector<diameter_index_t>& columns_to_reduce, entry_hash_map& pivot_column_index, index_t dim) {
        int cur_pos = next_simplex++;

        while (cur_pos < (int)simplices.size()) {
            diameter_index_t simplex = simplices[cur_pos];
            simplex_coboundary_enumerator cofacets(diameter_index_t(simplex), dim, *this);

			while (cofacets.has_next(false)) {
				auto cofacet = cofacets.next();
				if (get_diameter(cofacet) <= threshold) {

					outp.push_back(cofacet);

					if (pivot_column_index.find(get_index(cofacet)) == pivot_column_index.end())
						columns_to_reduce.push_back(cofacet);
				}
			}
			cur_pos = next_simplex++;
        }
	}

	void assemble_columns_to_reduce(std::vector<diameter_index_t>& simplices,
	                                std::vector<diameter_index_t>& columns_to_reduce,
	                                entry_hash_map& pivot_column_index, index_t dim, int num_workers=-1) {
        if (num_workers < 0)
            num_workers = std::thread::hardware_concurrency();

		--dim;
		columns_to_reduce.clear();
		std::vector<std::vector<diameter_index_t>> next_simplices(num_workers);
		std::vector<std::vector<diameter_index_t>> partitial_collumns(num_workers);

		std::vector<std::thread> ts;
		std::atomic<int> posit;
		posit = 0;

        for(int i = 0; i < num_workers; ++i)
            ts.emplace_back(&processor::parallel_assembling, this, std::ref(simplices), std::ref(posit),
                            std::ref(next_simplices[i]),  std::ref(partitial_collumns[i]),  std::ref(pivot_column_index), dim);

        int cofacets_num = 0;
        int collumns_num = 0;
        for (int i = 0; i < num_workers; ++i) {
            ts[i].join();
            cofacets_num += next_simplices[i].size();
            collumns_num += partitial_collumns[i].size();
        }

        simplices.clear();
        simplices.reserve(cofacets_num);
        columns_to_reduce.reserve(collumns_num);

        for (int i = 0; i < num_workers; ++i) {
            for (auto& x: next_simplices[i])
                simplices.push_back(x);
            for (auto& x: partitial_collumns[i])
                columns_to_reduce.push_back(x);
        }

		tbb::parallel_sort(columns_to_reduce.begin(), columns_to_reduce.end(),
		          compare_diameter_entries);
	}

	bool is_dominated(index_t s, index_t t, BitMatrix& working_graph) {
        auto common_vertices = working_graph[s] & working_graph[t];

        if (common_vertices.count() == 1)
            return true;

        for(int i = 0; i < (int)common_vertices.size(); ++i) {
            if (!common_vertices[i])
                continue;
            if ((common_vertices - working_graph[i]).count() == 1)
                return true;
        }
        return false;
	}

	void push_neighboorhood(BitMatrix& graph, BitMatrix& neighbourhood, index_t s, index_t t) {
        for (int row = 0; row < (int)neighbourhood.size(); ++row) {
            if (graph[row][s] & graph[row][t]) {
                neighbourhood[row][s] = neighbourhood[row][t] = 1;
                neighbourhood[s][row] = neighbourhood[t][row] = 1;
            }
        }
	}

    void delete_dominated_edges(std::vector<diameter_index_t>& edges, std::vector<bool>& dominated) {
        std::vector<index_t> vertices_of_edge(2);

        BitMatrix graph(n, boost::dynamic_bitset<>(n));
        BitMatrix backtrace_graph(n, boost::dynamic_bitset<>(n));
        BitMatrix neighbourhood(n, boost::dynamic_bitset<>(n));

        int i = 0;
        for (auto e : edges) {
			get_simplex_vertices(get_index(e), 1, n, vertices_of_edge.rbegin());
            graph[vertices_of_edge[0]][vertices_of_edge[1]] = graph[vertices_of_edge[1]][vertices_of_edge[0]] = 1;

			if (!is_dominated(vertices_of_edge[0], vertices_of_edge[1], graph))
                dominated[i] = false;
            else {
                ++i;
                continue;
            }

            backtrace_graph = graph;
            for (int row = 0; row < (int)neighbourhood.size(); ++row)
                neighbourhood[row].clear();

            push_neighboorhood(backtrace_graph, neighbourhood, vertices_of_edge[0], vertices_of_edge[1]);

            int j = i - 1;
            while (j >= 0) {
                get_simplex_vertices(get_index(edges[j]), 1, n, vertices_of_edge.rbegin());

                if (!dominated[j]) {
                    --j;
                    continue;
                }

                if (!neighbourhood[vertices_of_edge[0]][vertices_of_edge[1]]) {
                    backtrace_graph[vertices_of_edge[0]][vertices_of_edge[1]] = backtrace_graph[vertices_of_edge[1]][vertices_of_edge[0]] = 0;
                    --j;
                    continue;
                }

                if (!is_dominated(vertices_of_edge[0], vertices_of_edge[1], backtrace_graph)) {
                    dominated[j] = false;
                    push_neighboorhood(backtrace_graph, neighbourhood, vertices_of_edge[0], vertices_of_edge[1]);
                } else
                    backtrace_graph[vertices_of_edge[0]][vertices_of_edge[1]] = backtrace_graph[vertices_of_edge[1]][vertices_of_edge[0]] = 0;
                --j;
            }
            ++i;
		}
    }

	void compute_dim_0_pairs(std::vector<diameter_index_t>& meaningful_edges,
	                         std::vector<diameter_index_t>& columns_to_reduce,
	                         bool delete_dominated=true) {

		std::vector<diameter_index_t> edges = get_edges();
		std::sort(edges.rbegin(), edges.rend(), compare_diameter_entries);
		std::vector<index_t> vertices_of_edge(2);
        std::vector<bool> dominated(edges.size(), true);

        if (delete_dominated) {
            std::cerr << "Processing edge collapse... ";
	    auto  start = std::chrono::system_clock::now();
            delete_dominated_edges(edges, dominated);
            std::cerr << "done\n";
	    auto end_pos = std::chrono::system_clock::now();
	    std::chrono::duration<double> elapsed_seconds = end_pos - start;
	
	    std::cerr << "\n\nTime elapsed: " << elapsed_seconds.count() << " seconds\n\n";
        }

        std::cerr << "Searching for dim 0 pairs...";
        union_find dset(n);

		std::cout << "persistence intervals in dim 0:" << std::endl;

        int i = 0;
		for (auto e : edges) {
 			get_simplex_vertices(get_index(e), 1, n, vertices_of_edge.rbegin());

			if (delete_dominated && dominated[i]) {
                dist.rows[vertices_of_edge[1]][vertices_of_edge[0]] = threshold + 1.0;
                ++i;
                continue;
			}

            meaningful_edges.push_back(e);
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
				if (get_diameter(e) != 0)
					std::cout << " [0," << get_diameter(e) << ")" << std::endl;
				dset.link(u, v);
			} else
				columns_to_reduce.push_back(e);
            ++i;
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i)
                std::cout << " [0, )" << std::endl;
        std::cerr << "done\n";
	}

	template <typename Column> diameter_index_t pop_pivot(Column& column) {
		diameter_index_t pivot(0, -1);
		while (!column.empty()) {
			pivot = column.top();
			column.pop();
			if (column.empty() || get_index(column.top()) != get_index(pivot)) return pivot;
			column.pop();
		}
		return diameter_index_t(0, -1);
	}

	template <typename Column> diameter_index_t get_pivot(Column& column) {
		diameter_index_t result = pop_pivot(column);
		if (get_index(result) != -1) column.push(result);
		return result;
	}

	template <typename Column>
	diameter_index_t init_coboundary_and_get_pivot(const diameter_index_t simplex,
	                                               Column& working_coboundary, const index_t& dim,
	                                               entry_hash_map& pivot_column_index) {
		bool check_for_emergent_pair = true;
        thread_local static std::vector<diameter_index_t> cofacet_entries;
        cofacet_entries.clear();
		simplex_coboundary_enumerator cofacets(simplex, dim, *this);
		while (cofacets.has_next()) {
			diameter_index_t cofacet = cofacets.next();
			if (get_diameter(cofacet) <= threshold) {
				cofacet_entries.push_back(cofacet);
				if (check_for_emergent_pair && (get_diameter(simplex) == get_diameter(cofacet))) {
					if (pivot_column_index.find(get_index(cofacet)) == pivot_column_index.end())
						return cofacet;
					check_for_emergent_pair = false;
				}
			}
		}

		for (auto cofacet : cofacet_entries)
            working_coboundary.push(cofacet);

		return get_pivot(working_coboundary);
	}

	template <typename Column>
	void add_simplex_coboundary(const diameter_index_t simplex, const index_t& dim,
                                Column& working_reduction_column, Column& working_coboundary, bool add_diagonal=true) {

        if (add_diagonal)
            working_reduction_column.push(simplex);

		simplex_coboundary_enumerator cofacets(simplex, dim, *this);
		while (cofacets.has_next()) {
			diameter_index_t cofacet = cofacets.next();
			if (get_diameter(cofacet) <= threshold)
                working_coboundary.push(cofacet);
		}
	}

    template <typename Column>
	void add_coboundary(AsyncMatrix& reduction_matrix,
	                    const std::vector<diameter_index_t>& columns_to_reduce,
	                    const size_t index_column_to_add, const size_t& dim,
	                    PColumn loaded_column,
                        Column& working_reduction_column, Column& working_coboundary, bool add_diagonal=true) {

		diameter_index_t column_to_add(columns_to_reduce[index_column_to_add]);
		add_simplex_coboundary(column_to_add, dim, working_reduction_column, working_coboundary, add_diagonal);

        if (!loaded_column)
            return;

		for (diameter_index_t simplex : *loaded_column)
			add_simplex_coboundary(simplex, dim, working_reduction_column, working_coboundary);
	}

    template <class WorkingColumn>
    PColumn generate_column(WorkingColumn&& working_reduction_column) {
		if (working_reduction_column.empty())
			return nullptr;

		SparseColumn column;
		while (true) {
			diameter_index_t e = pop_pivot(working_reduction_column);
			if (get_index(e) == -1)
                break;

			column.push_back(e);
		}

		if (column.empty())
            return nullptr;
		return PColumn(new SparseColumn(std::move(column)));
	}

	void compute_pairs(const std::vector<diameter_index_t>& columns_to_reduce,
	                   entry_hash_map& pivot_column_index, const index_t dim, int num_workers=1) {

        std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
        AsyncMatrix reduction_matrix(columns_to_reduce.size());

        std::atomic<int> posit;
        posit = 0;
        std::vector<std::thread> ts;

        for (int i = 0; i < num_workers; ++i)
            ts.emplace_back(&processor::parallel_reduction, this, std::ref(reduction_matrix), std::ref(columns_to_reduce), std::ref(pivot_column_index),
                            std::ref(posit), dim, i);

        for (int i = 0; i < num_workers; ++i)
            ts[i].join();

        for (auto& persistance_pair: pivot_column_index)
            if (persistance_pair.second.second > get_diameter(columns_to_reduce[persistance_pair.second.first]))
                std::cout << " [" << get_diameter(columns_to_reduce[persistance_pair.second.first]) << "," << persistance_pair.second.second << ")\n";
    }

    void parallel_reduction(AsyncMatrix& reduction_matrix, const std::vector<diameter_index_t>& columns_to_reduce,
	                   entry_hash_map& pivot_column_index, std::atomic<int>& next_column,
	                   const index_t dim, int id) {

        int index_column_to_reduce = next_column++;
        int nxt;
        bool first;

        while (index_column_to_reduce < (int)reduction_matrix.size()) {
            first = true;
            do {
                nxt = index_column_to_reduce;
                index_column_to_reduce = one_reduction_step(reduction_matrix, columns_to_reduce, pivot_column_index, nxt, dim, first);
                first = false;
            } while (nxt != index_column_to_reduce);

            index_column_to_reduce = next_column++;
        }
    }

	index_t one_reduction_step(AsyncMatrix& reduction_matrix, const std::vector<diameter_index_t>& columns_to_reduce,
	                   entry_hash_map& pivot_column_index, int index_column_to_reduce, const index_t dim, bool first) {

        diameter_index_t column_to_reduce(columns_to_reduce[index_column_to_reduce]);

        std::priority_queue<diameter_index_t, std::vector<diameter_index_t>, greater_diameter_or_smaller_index>
                            working_reduction_column, working_coboundary;

        diameter_index_t pivot;

        if (first)
            pivot = init_coboundary_and_get_pivot(column_to_reduce, working_coboundary, dim, pivot_column_index);
        else {
            PColumn initial_column = std::atomic_load(&reduction_matrix[index_column_to_reduce]);
            add_coboundary(reduction_matrix, columns_to_reduce, index_column_to_reduce, dim, initial_column,
                            working_reduction_column, working_coboundary, false);
            pivot = get_pivot(working_coboundary);
        }

        while (true) {
            if (get_index(pivot) != -1) {
                auto pair = pivot_column_index.find(get_index(pivot));
                if (pair != pivot_column_index.end()) {

                    index_t index_column_to_add = pair->second.first;;
                    index_t old_index_column_to_add;
                    PColumn column_to_add = nullptr;

                    do {
                        old_index_column_to_add = index_column_to_add;
                        column_to_add = std::atomic_load(&reduction_matrix[index_column_to_add]);

                        index_column_to_add = pair->second.first;
                    } while (old_index_column_to_add != index_column_to_add);

                    if (index_column_to_add < index_column_to_reduce) {
                        add_coboundary(reduction_matrix, columns_to_reduce, index_column_to_add, dim, column_to_add,
                                        working_reduction_column, working_coboundary);
                        pivot = get_pivot(working_coboundary);

                    } else {
                        PColumn new_column = generate_column(std::move(working_reduction_column));
                        std::atomic_store(&reduction_matrix[index_column_to_reduce], new_column);

                        if (pair->second.first.compare_and_swap(index_column_to_reduce, index_column_to_add)) {
                            return index_column_to_add;
                        } else {
                            continue;
                        }
                    }
                } else {
                    PColumn new_column = generate_column(std::move(working_reduction_column));
                    std::atomic_store(&reduction_matrix[index_column_to_reduce], new_column);

                    value_t death = get_diameter(pivot);

                    auto insertion_result = pivot_column_index.insert({get_index(pivot), {index_column_to_reduce, death}});

                    if (!insertion_result.second) {
                        continue;
                    }
                    break;
                }
            } else{
                break;
            }
        }
        return index_column_to_reduce;
    }

	std::vector<diameter_index_t> get_edges();

	void compute_barcodes(bool use_collapse, int num_workers = -1) {
        if (num_workers < 0)
            num_workers = std::thread::hardware_concurrency();

		std::vector<diameter_index_t> simplices, columns_to_reduce;

		compute_dim_0_pairs(simplices, columns_to_reduce, use_collapse);

		for (index_t dim = 1; dim <= dim_max; ++dim) {
			entry_hash_map pivot_column_index;

			std::cerr << "Processing dim " << char(dim +'0') << " pairs... ";
			compute_pairs(columns_to_reduce, pivot_column_index, dim, num_workers);
			std::cerr << "done";

			if (dim < dim_max) {
                std::cerr << ". Assembling columns for the next step... ";
				assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
				                           dim + 1, num_workers);
                std::cerr << "done";
            }
            std::cerr <<"\n";
		}
	}
};

std::vector<diameter_index_t> processor::get_edges() {
	std::vector<diameter_index_t> edges;
	std::vector<index_t> vertices(2);

	for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
		get_simplex_vertices(index, 1, dist.size(), vertices.rbegin());

		value_t length = dist(vertices[0], vertices[1]);
		if (length <= threshold && length >= 0)
            edges.push_back({length, index});
	}
	return edges;
}

compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;

	std::string line;
	value_t value;
	for (int i = 0; std::getline(input_stream, line); ++i) {
		std::istringstream s(line);
		for (int j = 0; j < i && s >> value; ++j) {
			distances.push_back(value);
			s.ignore();
		}
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

value_t euclidean_distance(std::vector<value_t>& point_1, std::vector<value_t>& point_2) {
    value_t ans = 0.0;

    if (point_1.size() != point_2.size()) {
        std::cerr << "Invalid input data format. Please check dimension of data points\n";
        exit(-1);
    }

    for (int i = 0; i < (int)point_1.size(); ++i)
        ans += (point_1[i] - point_2[i]) * (point_1[i] - point_2[i]);

    return sqrt(ans);
}

compressed_lower_distance_matrix read_point_cloud(std::istream& input_stream) {
    std::vector<std::vector<value_t>> points_cloud;
	std::vector<value_t> distances;

	std::string line;
	value_t value;
	for (int i = 0; std::getline(input_stream, line); ++i) {
		std::istringstream s(line);

		points_cloud.push_back(std::vector<value_t>());
		for(; s >> value; ) {
            points_cloud[i].push_back(value);
            s.ignore();
		}
    }

    for (int i = 0; i < (int)points_cloud.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            distances.push_back(euclidean_distance(points_cloud[i], points_cloud[j]));
        }
    }

	return compressed_lower_distance_matrix(std::move(distances));
}

void exit_with_info(int exit_code) {
	std::cerr
	    << "Usage: "
	    << "PersistenceCalculator "
	    << "[input file name] [options]\n\n"
	    << "Options:\n\n"
	    << "  --help           print this screen\n"
        << "  --output         output file name (DEFAULT: output.txt )\n"
	    << "  --format         use the specified file format for the input. Options are:\n"
	    << "                     matrix         (distance matrix, DEFAULT)\n"
	    << "                     cloud          (point cloud in Euclidean space)\n"
	    << "  --dim <k>        compute persistent homology up to dimension k (DEFAULT: 2)\n"
	    << "  --threshold <t>  compute Rips complexes up to diameter t\n"
	    << "  --prep           force to use preprocessing by eliminating dominated edges\n"
	    << "  --noprep         force to NOT use preprocessing\n"
	    << "  --threads <n>    use n threads for multithreading (DEFAULT: per CPU available)\n";
	exit(exit_code);
}

int main(int argc, char** argv) {

	const char* output = "output.txt";
	index_t dim_max = 2;
	value_t threshold = std::numeric_limits<value_t>::max();

	int input_type = 0;
	int preprocess_type = 0;
	int num_workers = -1;

	if (std::string(argv[1]) == "--help")
        exit_with_info(0);

    const char* filename = argv[1];
	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	for (index_t i = 2; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--dim") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			dim_max = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size())
                exit_with_info(-1);
		} else if (arg == "--threshold") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			threshold = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size())
                exit_with_info(-1);
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter.rfind("cloud", 0) == 0)
			    input_type = 1;
			else if (parameter.rfind("matrix", 0) == 0)
			    input_type = 0;
			else
			    exit_with_info(-1);
		} else if (arg == "--threads") {
            		std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			num_workers = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size())
                exit_with_info(-1);
		} else if (arg == "--output") {
            output = argv[++i];
            if (!output)
                exit_with_info(-1);
		} else if (arg == "--prep") {
            preprocess_type = 1;
		} else if (arg == "--noprep") {
            preprocess_type = -1;
		} else
		{
            exit_with_info(-1);
		}
	}

    compressed_lower_distance_matrix dist;
    if (input_type == 0)
        dist = read_distance_matrix(filename ? file_stream : std::cin);
    else
        dist = read_point_cloud(filename ? file_stream : std::cin);

    bool use_edge_collapse = true;

    if (preprocess_type == -1)
        use_edge_collapse = false;
    if (preprocess_type == 0)
        use_edge_collapse = (dim_max >= 3);

    freopen(output, "w", stdout);

    auto start = std::chrono::system_clock::now();

    value_t min = std::numeric_limits<value_t>::infinity(), max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
    int num_edges = 0;

    if (threshold == std::numeric_limits<value_t>::max()) {
        value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
        for (size_t i = 0; i < dist.size(); ++i) {
            value_t r_i = -std::numeric_limits<value_t>::infinity();

            for (size_t j = 0; j < dist.size(); ++j)
                r_i = std::max(r_i, dist(i, j));

            enclosing_radius = std::min(enclosing_radius, r_i);
        }
        threshold = enclosing_radius;
    }

    for (auto d : dist.distances) {
        min = std::min(min, d);
        max = std::max(max, d);
        max_finite = (d != std::numeric_limits<value_t>::infinity() ? std::max(max, d) : max_finite);
        if (d <= threshold)
            ++num_edges;
    }

    std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;
    processor(std::move(dist), dim_max, threshold).compute_barcodes(use_edge_collapse, num_workers);

    auto end_pos = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_pos - start;

    std::cerr << "\n\nTime elapsed: " << elapsed_seconds.count() << " seconds\n";

    return 0;
}
