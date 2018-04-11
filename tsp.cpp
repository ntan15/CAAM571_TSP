/*
 * TSP.cpp
 *
 * Nigel, Nick, Shengchao
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <time.h>
#include <limits>
#include <algorithm>
#include <utility>
#include <tuple>
#include <map>
#include <unordered_map>
#include <functional>
#include "gurobi_c++.h"

//#define UNORDERED

// Ignore any hash stuff. 
/*
struct pairhash {
    std::size_t operator()(const std::pair<unsigned int,unsigned int>& x) const {
        unsigned int h1 = ((x.first >> 16) ^ x.first) * 0x45d9f3b;
        h1 = ((h1 >> 16) ^ h1) * 0x45d9f3b;
        h1 = (h1 >> 16) ^ h1;
        unsigned int h2 = ((x.second >> 16) ^ x.second) * 0x45d9f3b;
        h2 = ((h2 >> 16) ^ h2) * 0x45d9f3b;
        h2 = (h2 >> 16) ^ h2;
        return h1 ^ h2;
    }
};
*/
struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

class abc {
  std::unordered_map<std::pair<int,int>, int, pairhash> rules;
};

using namespace std;

unsigned int hash_test(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

/*
class Edges {
    private:
        vector< array<int,3> > e;
        vector<int> start_idx;
    public:
        Edges(vector< array<int,3> >&& vec, vector<int>&& idx);
        int get_weight(int a, int b);
};
*/
class Partition {
    private:
        vector< vector<int> > v;
    public:
        Partition(vector< vector<int> > && vec);
        int contains(int val);
//        int distance(Partition& B, Edges& e);
        int num_groups();
        int num_nodes();
        void add_group(vector<int>&& group);
        void remove_group(int idx);
        vector<int> get_group(int idx);
        void add_to_group(int group, vector<int>& val);
        void print();
};
int contains_val(vector< vector<int> >& v, int val);
//int min_cut_phase(Partition& A, Partition& B, Edges e, int a);
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, vector<array<int,3> >& e, int a);
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e, int a);
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, unordered_map<pair<int,int>,int>& e, int a);
//int min_cut(vector<int>& v, Edges& e, int a);
tuple<int,vector<pair<int,int> > > min_cut(vector<int>& v, vector<array<int,3> >& e, int a);
tuple<int,vector<pair<int,int> > > min_cut(vector<int>& v, map<pair<int,int>,int>& e, int a);
tuple<int,vector<pair<int,int> > > min_cut(vector<int>& v, unordered_map<pair<int,int>,int>& e, int a);
int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, vector<array<int,3> >& e);
int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e);
int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, unordered_map<pair<int,int>,int>& e);
int contains_val(vector< vector<int> >& v, int val);

/*
Edges::Edges(vector< array<int,3> >&& vec, vector<int>&& idx)
{
    e = vec;
    start_idx = idx;
}
int Edges::get_weight(int a, int b)
{
// Naive
    for(int i=0; i<e.size(); i++) {
        if( ((e[i][0]==a) && (e[i][1]==b)) || ((e[i][0]==b) && (e[i][1]==a)))
        {
            return e[i][2];
        }
    }
    return 0;
}
*/

Partition::Partition(vector< vector<int> > && vec)
{
    v = vec;
}
int Partition::contains(int val)
{
    for(int i=0; i<v.size(); i++)
    {
        if(std::find(v[i].begin(),v[i].end(),val) != v[i].end())
        {
            return i;
        }
    }
    return -1;       
}
/*
int Partition::distance(Partition& B, Edges& e)
{
    int dist = 0;
    for(int i=0; i<v.size(); i++)
    {
        for(int j=0; j<v[i].size(); j++)
        {
            for(int k=0; k<B.v.size(); k++)
            {
                for(int l=0; l<B.v[i].size(); l++)
                {
                    dist += e.get_weight(v[i][j], B.v[i][j]);
                }
            }
        }
    }
    return dist;
}
*/
int Partition::num_groups() {
    return v.size();
}
int Partition::num_nodes() {
    int size = 0;
    for(int i=0; i<v.size(); i++) {
        size += v[i].size();
    }
    return size;
}
void Partition::add_group(vector<int>&& group) {
    v.push_back(group);
}
void Partition::remove_group(int idx) {
    v.erase(v.begin()+idx);
}
vector<int> Partition::get_group(int idx) {
    return v[idx];
}
void Partition::add_to_group(int group, vector<int>& val) {
    v[group].insert(v[group].end(), val.begin(), val.end());
}
void Partition::print() {
    for(int j=0; j<v.size(); j++) {
        cout << "(";
        for(int k=0; k<v[j].size(); k++) {
            cout << v[j][k] << " ";
        }
        cout << "),";
    }
    cout << endl;
}
/*
int min_cut_phase(Partition& A, Partition& B, Edges e, int a) {
    while(B.num_groups() != 2) {
        int tight_idx = most_tight(A,B,e);
        A.add_group(std::move(B.get_group(tight_idx)));
        B.remove_group(tight_idx);
    }
    B.add_to_group(0,B.get_group(1));
    V.remove_group(1);
}
*/
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, vector<array<int,3> >& e, int a)
{
    while(B.size() > 2)
    {
        int tight_idx = most_tight(A,B,e);
        A.push_back(B[tight_idx]);
        B.erase(B.begin()+tight_idx);
    }
    B[0].insert(B[0].end(), B[1].begin(), B[1].end());
    B.erase(B.begin()+1);
// Generate cut
    vector<pair<int,int> > edges;
    int weight = 0;
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            for(int k=0; k<B.size(); k++) {
                for(int l=0; l<B[k].size(); l++) {
                    for(int p=0; p<e.size(); p++) {
                        if( ((e[p][0]==A[i][j]) && (e[p][1]==B[k][l])) || ((e[p][0]==B[k][l]) && (e[p][1]==A[i][j])))
                        {
                            weight += e[p][2];
                            edges.push_back(make_pair(e[p][0],e[p][1]));
                        }
                    }
                }
            }
        }
    }
    return make_tuple(weight,edges);
}
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e, int a)
{
    while(B.size() > 2)
    {
        int tight_idx = most_tight(A,B,e);
        A.push_back(B[tight_idx]);
        B.erase(B.begin()+tight_idx);
    }
    B[0].insert(B[0].end(), B[1].begin(), B[1].end());
    B.erase(B.begin()+1);
// Generate cut
    vector<pair<int,int> > edges;
    int weight = 0;
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            for(int k=0; k<B.size(); k++) {
                for(int l=0; l<B[k].size(); l++) {
                    auto iter = e.find(pair<int,int>(A[i][j],B[k][l]));
                    if(iter != e.end()) {
                        weight += iter->second;
                        edges.push_back(make_pair(A[i][j],B[k][l]));
                    } else {
                        iter = e.find(pair<int,int>(B[k][l],A[i][j]));
                        if(iter != e.end()) {
                            weight += iter-> second;
                            edges.push_back(make_pair(B[k][l],A[i][j]));
                        }
                    }
                }
            }
        }
    }
    return make_tuple(weight,edges);
}
#ifdef UNORDERED
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, unordered_map<pair<int,int>,int>& e, int a)
{
    while(B.size() > 2)
    {
        int tight_idx = most_tight(A,B,e);
        A.push_back(B[tight_idx]);
        B.erase(B.begin()+tight_idx);
    }
    B[0].insert(B[0].end(), B[1].begin(), B[1].end());
    B.erase(B.begin()+1);
// Generate cut
    vector<pair<int,int> > edges;
    int weight = 0;
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            for(int k=0; k<B.size(); k++) {
                for(int l=0; l<B[k].size(); l++) {
                    auto iter = e.find(pair<int,int>(A[i][j],B[k][l]));
                    if(iter != e.end()) {
                        weight += iter->second;
                        edges.push_back(make_pair(A[i][j],B[k][l]));
                    } else {
                        iter = e.find(pair<int,int>(B[k][l],A[i][j]));
                        if(iter != e.end()) {
                            weight += iter-> second;
                            edges.push_back(make_pair(B[k][l],A[i][j]));
                        }
                    }
                }
            }
        }
    }
    return make_tuple(weight,edges);
}
#endif
/*
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e, int a)
{
    while(B.size() > 2)
    {
        int tight_idx = most_tight(A,B,e);
cout << "Tightest node in B: B[" << tight_idx << "] = ";
for(int i=0; i<B[tight_idx].size(); i++) {
cout << B[tight_idx][i] << " ";
}
cout << endl;
        A.push_back(B[tight_idx]);
        B.erase(B.begin()+tight_idx);
cout << "Merged Vertices" << endl;
    }
cout << "Creating Final partitions\n";
    B[0].insert(B[0].end(), B[1].begin(), B[1].end());
    B.erase(B.begin()+1);
cout << "Created Final partitions\n";
// Generate cut
cout << "Compute distance between partitions\n";
    vector<pair<int,int> > edges;
    int weight = 0;
    pair<int,int> e_idx;
    map<pair<int,int>,int>::iterator iter;
    
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            for(int k=0; k<B.size(); k++) {
                for(int l=0; l<B[k].size(); l++) {
                    e_idx = pair<int,int>(A[i][j],B[k][l]);
                    iter = e.find(e_idx);
                    if(iter != e.end()) {
                        weight += *iter;
                        edges.push_back(e_idx);
                    } else {
                        e_idx = pair<int,int>(B[k][l],A[i][j]);
                        iter = e.find(e_idx);
                        if(iter != e_end()) {
                            weight += *iter;
                            edges.push_back(e_idx);
                        }
                    }
                }
            }
        }
    }
cout << "Computed distance between partitions\n";
    return make_tuple(weight,edges);
}
*/
/*
int min_cut(vector<int>& v, Edges& e, int a) {
    int min_cut = numeric_limits<int>::max();
    int current_cut;
}
*/
tuple<int,vector<pair<int,int> > > min_cut(vector<int>& v, vector<array<int,3> >& e, int a)
{
    int min_cut = numeric_limits<int>::max();
    int current_weight;
    vector<pair<int,int> > current_cut;
    vector< vector<int> > A;
    A.push_back(vector<int>(1,a));
    vector< vector<int> > B;
    for(int i=0; i<v.size(); i++)
    {
        if(v[i] != a)
            B.push_back(vector<int>(1,v[i]));
    }
    for(int i=0; i<v.size()-1; i++)
    {
        auto cut  = min_cut_phase(A,B,e,a);
        current_weight = get<0>(cut);
        if( (current_weight < min_cut) && (current_weight > 0) )
        {
            min_cut = current_weight;
            current_cut = get<1>(cut);
        }
        B.insert(B.end(),A.begin()+1, A.end());
        A.erase(A.begin()+1,A.end());
    }
    return make_tuple(min_cut,current_cut);
}
tuple<int,vector<pair<int,int> > > min_cut(vector<int>& v, map<pair<int,int>,int>& e, int a)
{
    int min_cut = numeric_limits<int>::max();
    int current_weight;
    vector<pair<int,int> > current_cut;
    vector< vector<int> > A;
    A.push_back(vector<int>(1,a));
    vector< vector<int> > B;
    for(int i=0; i<v.size(); i++)
    {
        if(v[i] != a)
            B.push_back(vector<int>(1,v[i]));
    }
    for(int i=0; i<v.size()-1; i++)
    {
        auto cut  = min_cut_phase(A,B,e,a);
        current_weight = get<0>(cut);
        if( (current_weight < min_cut) && (current_weight > 0) )
        {
            min_cut = current_weight;
            current_cut = get<1>(cut);
        }
        B.insert(B.end(),A.begin()+1, A.end());
        A.erase(A.begin()+1,A.end());
    }
    return make_tuple(min_cut,current_cut);
}
#ifdef UNORDERED
tuple<int,vector<pair<int,int> > > min_cut(vector<int>& v, unordered_map<pair<int,int>,int>& e, int a)
{
    int min_cut = numeric_limits<int>::max();
    int current_weight;
    vector<pair<int,int> > current_cut;
    vector< vector<int> > A;
    A.push_back(vector<int>(1,a));
    vector< vector<int> > B;
    for(int i=0; i<v.size(); i++)
    {
        if(v[i] != a)
            B.push_back(vector<int>(1,v[i]));
    }
    for(int i=0; i<v.size()-1; i++)
    {
        auto cut  = min_cut_phase(A,B,e,a);
        current_weight = get<0>(cut);
        if( (current_weight < min_cut) && (current_weight > 0) )
        {
            min_cut = current_weight;
            current_cut = get<1>(cut);
        }
        B.insert(B.end(),A.begin()+1, A.end());
        A.erase(A.begin()+1,A.end());
    }
    return make_tuple(min_cut,current_cut);
}
#endif

int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, vector<array<int,3> >& e)
{
    int idx = -1;
    int start_in_A, end_in_A, start_in_B, end_in_B;
    vector<int> sum_weights = vector<int>(B.size(),0);
    for(int j=0; j<e.size(); j++)
    {
        start_in_A = contains_val(A,e[j][0]);
        end_in_B = contains_val(B,e[j][1]);
        end_in_A = contains_val(A,e[j][1]);
        start_in_B = contains_val(B,e[j][0]);
        if((start_in_A!=-1) && (end_in_B!=-1))
        {
            sum_weights[end_in_B] += e[j][2];
        }
        if((start_in_B!=-1) && (end_in_A!=-1))
        {
            sum_weights[start_in_B] += e[j][2];
        }
    }
    idx = max_element(sum_weights.begin(), sum_weights.end())-sum_weights.begin();
    return idx;
}
int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e)
{
    int idx = -1;
    int start_in_A, end_in_A, start_in_B, end_in_B;
    vector<int> sum_weights = vector<int>(B.size(),0);
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            for(int k=0; k<B.size(); k++) {
                for(int l=0; l<B[k].size(); l++) {
                    auto iter = e.find(pair<int,int>(A[i][j],B[k][l]));
                    if(iter != e.end()) {
                        sum_weights[k] += iter->second;
                        continue;
                    }
                    iter = e.find(pair<int,int>(B[k][l],A[i][j]));
                    if(iter != e.end()) {
                        sum_weights[k] += iter->second;
                    }
                }
            }
        }
    }
    idx = max_element(sum_weights.begin(), sum_weights.end())-sum_weights.begin();
    return idx;
}
#ifdef UNORDERED
int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, unordered_map<pair<int,int>,int>& e)
{
    int idx = -1;
    int start_in_A, end_in_A, start_in_B, end_in_B;
    vector<int> sum_weights = vector<int>(B.size(),0);
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            for(int k=0; k<B.size(); k++) {
                for(int l=0; l<B[k].size(); l++) {
                    auto iter = e.find(pair<int,int>(A[i][j],B[k][l]));
                    if(iter != e.end()) {
                        sum_weights[k] += iter->second;
                        continue;
                    }
                    iter = e.find(pair<int,int>(B[k][l],A[i][j]));
                    if(iter != e.end()) {
                        sum_weights[k] += iter->second;
                    }
                }
            }
        }
    }
    idx = max_element(sum_weights.begin(), sum_weights.end())-sum_weights.begin();
    return idx;
}
#endif

int contains_val(vector< vector<int> >& v, int val)
{
    for(int i=0; i<v.size(); i++)
    {
        if(std::find(v[i].begin(),v[i].end(),val) != v[i].end())
        {
            return i;
        }
    }
    return -1;       
}


int main(int argc, char* argv[]) {
    ifstream f_in;
    ofstream f_out;
    int n; // # of nodes
    int m; // # of edges
    vector< array<int,3> > edges;
    map<pair<int,int>,int> edge_map;
pairhash test_h;
    unordered_map<pair<int,int>,int,pairhash> unordered_edge_map;
    vector<int> nodes;

    if(argc < 2) {
        cerr << "Usage: tsp input_file.txt" << endl;
        return 0;
    }
    f_in.open(argv[1]);
    if(!f_in.is_open()) {
        cerr << "Could not open tsp file!" << endl;
        return -1;
    }
// Read node/edge sizes
    f_in >> n >> m;
// Read edges
    for(int i=0; i<m; i++) {
        int start, end, w;
        f_in >> start >> end >> w;
        array<int,3> arr;
        arr[0] = start;
        arr[1] = end;
        arr[2] = w;
        edges.push_back(arr);
        edge_map[pair<int,int>(start,end)] = w;
        unordered_edge_map[pair<int,int>(start,end)] = w;
    }
// Set Nodes
    for(int i=0; i<n; i++) {
        nodes.push_back(i);
    }
abc test_hash_struct;
/*
    try {
        // Declare the GRB environment
        GRBEnv env = GRBEnv();
        // set the model within the GRB environment
        GRBModel model= GRBModel(env);
        vector<GRBVar> x;
        for(int i=0; i<m; i++) {
            x.push_back(model.addVar(0,1,0,GRB_CONTINUOUS));
        }
        model.update();

    } catch(GRBException e) {
        cout << "Gurobi Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during optimization" << endl;
    }
*/
clock_t timer = clock();
tuple<int,vector<pair<int,int> > > test_map = min_cut(nodes,edge_map,1);
timer = clock() - timer;
cout.precision(numeric_limits<double>::max_digits10);
cout << "Min Cut Algorithm Map Time: " << scientific << ((double)timer)/(double)CLOCKS_PER_SEC << endl;
cout << "Min Cut Algorithm Map Clock Tics: " << timer << endl;
int weight = get<0>(test_map);
vector<pair<int,int> > min_edges = get<1>(test_map);
cout << "Cost of Min cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

timer = clock();
tuple<int,vector<pair<int,int> > > test = min_cut(nodes, edges, 1);
timer = clock() - timer;
cout.precision(numeric_limits<double>::max_digits10);
cout << "Min Cut Algorithm Time: " << scientific << ((double)timer)/(double)CLOCKS_PER_SEC << endl;
cout << "Min Cut Algorithm Clock Tics: " << timer << endl;
weight = get<0>(test);
min_edges = get<1>(test);
cout << "Cost of Min cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}
/*
timer = clock();
test = min_cut(nodes, unordered_edge_map, 1);
timer = clock() - timer;
cout.precision(numeric_limits<double>::max_digits10);
cout << "Min Cut Algorithm Time: " << scientific << ((double)timer)/(double)CLOCKS_PER_SEC << endl;
cout << "Min Cut Algorithm Clock Tics: " << timer << endl;
weight = get<0>(test);
min_edges = get<1>(test);
cout << "Cost of Min cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}
*/
cout << "hash(1): " << hash_test(1) << endl;
cout << "hash(2): " << hash_test(2) << endl;
cout << "hash(3): " << hash_test(3) << endl;
cout << "hash(4): " << hash_test(4) << endl;
cout << "hash(6): " << hash_test(6) << endl;
cout << "hash_test(2) ^ hash_test(4): " << (hash_test(2) ^ hash_test(4)) << endl;
cout << "hash_test(4) ^ hash_test(2): " << (hash_test(4) ^ hash_test(2)) << endl;
cout << "hash_test(3) ^ hash_test(3): " << (hash_test(3) ^ hash_test(3)) << endl;
cout << "hash_test(2) ^ hash_test(3): " << (hash_test(2) ^ hash_test(3)) << endl;
cout << "hash_test(6) ^ hash_test(1): " << (hash_test(6) ^ hash_test(1)) << endl;
cout << "Unordered Map: " << (*(unordered_edge_map.find(pair<int,int>(0,1)))).second << endl;

// Output Results
    f_out.open("test_output.out");
    for(int i=0; i<m; i++) {
        f_out << edges[i][0] << " " << edges[i][1] << " " << edges[i][2] << " //for edge(" << i << ") of tour" << endl;
    }
    f_out << "The cost of the best tour is: (the cost the best tour" << endl;
    f_out << "Cost of Min cut: " << weight << endl;
    for(int i=0; i<min_edges.size(); i++) {
        f_out << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
    }
    f_out.close();
    return 1;
}

