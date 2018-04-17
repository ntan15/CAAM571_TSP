/*
 * TSP.cpp
 *
 * Nigel, Nick, Shengchao
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include <limits>
#include <algorithm>
#include <utility>
#include <tuple>
#include <map>
#include <stack>
#include <array>
#include <functional>
#include <chrono>
#include "gurobi_c++.h"

#define STACK

using namespace std;

unsigned int hash_test(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

struct pair_hash {
    template <class T1>
    size_t operator() (const pair<T1,T1>& p) const {
        auto h1 = hash_test(p.first);
        auto h2 = hash_test(p.second);
        return h1 ^ h2;
    }
};

bool same(double* x, double* y, int m);
tuple<double,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,double>& e, int a);
tuple<double,vector<pair<int,int> > > min_cut(vector<int>& v, map<pair<int,int>,double>& e, int a);
int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,double>& e);
bool is_binary(map<pair<int,int>,GRBVar>& x, int m);
double tsp(GRBModel& model, map<pair<int,int>,GRBVar>& x, pair<GRBLinExpr,GRBLinExpr>& constr, vector<double>& best_x, double& best_cost, int n, int m, vector<int>& nodes);
double tsp(GRBModel& model, map<pair<int,int>,GRBVar>& x, pair<GRBLinExpr,GRBLinExpr>& constr, vector<double>& best_x, double& best_cost, int n, int m, vector<int>& nodes, vector<array<double,3> >& pseudo_cost, int prior_branch, double prior_cost, double prior_x);
double nearest_neighbor(vector<int>& V, map<pair<int,int>, double>& E);

bool same(double* x, double* y, int m) {
    for(int i=0; i<m; i++) {
        if(abs(x[i] - y[i]) > pow(10,-7)) {
            return false;
        }
    }
    return true;
}

double nearest_neighbor(vector<int>& V, map<pair<int,int>, double>& E, int a) {
    int n = V.size();
    vector<bool> visited = vector<bool>(V.size(), false);
    int start = a;
    int origin = start;
    int end = start;
    double dist = 0;
    for(int i=0; i<n-1; i++) {
        visited[start] = true;
        double nearest = numeric_limits<double>::max();
        for(int j=0; j<n; j++) {
            if(!visited[j]) {
                auto edge = E.find(pair<int,int>(start,j));
                if( (edge != E.end()) && (edge->second < nearest) ) {
                    nearest = edge->second;
                    end = j;
                }
                edge = E.find(pair<int,int>(j,start));
                if( (edge != E.end()) && (edge->second < nearest) ) {
                    nearest = edge->second;
                    end = j;
                }
            }
        }
        start = end;
        dist += nearest;
    }
    auto edge = E.find(pair<int,int>(start,origin));
    if(edge != E.end()) {
        dist += edge->second;
    } else {
        edge = E.find(pair<int,int>(origin,start));
        dist += edge->second;
    }
    return dist;
}

tuple<double,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,double>& e, int a)
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
    double weight = 0;
    int edge_count = 0;
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<A[i].size(); j++) {
            for(int k=0; k<B.size(); k++) {
                for(int l=0; l<B[k].size(); l++) {
                    auto iter = e.find(pair<int,int>(A[i][j],B[k][l]));
                    if(iter != e.end()) {
                        edge_count++;
                        weight += iter->second;
                        edges.push_back(make_pair(A[i][j],B[k][l]));
                    } else {
                        iter = e.find(pair<int,int>(B[k][l],A[i][j]));
                        if(iter != e.end()) {
                            edge_count++;
                            weight += iter->second;
                            edges.push_back(make_pair(B[k][l],A[i][j]));
                        }
                    }
                }
            }
        }
    }
if(edge_count == 0)
weight = -1;
    return make_tuple(weight,edges);
}

tuple<double,vector<pair<int,int> > > min_cut(vector<int>& v, map<pair<int,int>,double>& e, int a)
{
    double min_cut = numeric_limits<double>::max();
    double current_weight;
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
        if( (current_weight < min_cut) && (current_weight!=-1))
        {
            min_cut = current_weight;
            current_cut = get<1>(cut);
        }
        B.insert(B.end(),A.begin()+1, A.end());
        A.erase(A.begin()+1,A.end());
    }
    return make_tuple(min_cut,current_cut);
}

int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,double>& e)
{
    int idx = -1;
    vector<double> sum_weights = vector<double>(B.size(),0);
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

bool is_binary(map<pair<int,int>,GRBVar>& x) {
    for(auto it=x.begin(); it!=x.end(); ++it) {
        double xval = (it->second).get(GRB_DoubleAttr_X);
        if( abs(xval - round(xval)) > pow(10,-10) ) {
            return false;
        }
    }
    return true;
}

double tsp(GRBModel& model, map<pair<int,int>,GRBVar>& x, pair<GRBLinExpr,GRBLinExpr>& constr, 
           vector<double>& best_x, double& best_cost, int n, int m, vector<int>& nodes, vector<array<double,3> >& pseudo_cost, int prior_branch, double prior_cost, double prior_x) {
    auto constraint = model.addConstr(constr.first == constr.second);
    model.optimize();
    vector<double> x_val = vector<double>(m,0);
    map<pair<int,int>,double> lp_edge_map;
    
    if(model.get(GRB_IntAttr_Status) == 3) {
        model.remove(constraint);
        return best_cost;
    }

    double C = model.get(GRB_DoubleAttr_ObjVal);


    if(C >= best_cost) {
        model.remove(constraint);
        return best_cost;
    } else {
        for(int i=0; i<m; i++) {
            x_val[i] = (next(x.begin(),i)->second).get(GRB_DoubleAttr_X);
            lp_edge_map[next(x.begin(),i)->first] = x_val[i];
        }
        unsigned int frac = -1;
        double frac_diff = 1.0;
        double down_est = 0;
        double up_est = 0;
        double best_gain = 0;
        for(int i=0; i<m; i++) {
            if(abs(x_val[i] - 0.5) < frac_diff) {
                frac_diff = abs(x_val[i] - 0.5);
                frac = i;
            }
        }

        auto mincut = min_cut(nodes,lp_edge_map,rand()%n);
        double weight = get<0>(mincut);
        if(abs(round(weight)-weight) < pow(10,-10)) {
            weight = round(weight);
        }

        while( weight < 2 ) {
        // Subtour
            GRBLinExpr st_expr = 0;
            auto mincut_edges = get<1>(mincut);
            mincut_edges = get<1>(mincut);
            for(int i=0; i<mincut_edges.size(); i++) {
                st_expr += x[mincut_edges[i]];
            }
            model.addConstr(st_expr >= 2);
            model.optimize();
            if(model.get(GRB_IntAttr_Status) == 3) {
                model.remove(constraint);
                return best_cost;
            }
            C = model.get(GRB_DoubleAttr_ObjVal);
            if(C >= best_cost) {
                model.remove(constraint);
                return best_cost;
            }
            for(int i=0; i<m; i++) {
                x_val[i] = (next(x.begin(),i)->second).get(GRB_DoubleAttr_X);
                lp_edge_map[next(x.begin(),i)->first] = x_val[i];
            }
            frac = -1;
            frac_diff = 1.0;
            for(int i=0; i<m; i++) {
                if( (abs(round(x_val[i]) - x_val[i]) > pow(10,-10)) ) {
                    if(abs(x_val[i] - 0.5) < frac_diff) {
                        frac_diff = abs(x_val[i] - 0.5);
                        frac = i;
                    }
                }
            }
            mincut = min_cut(nodes,lp_edge_map,rand()%n);
            weight = get<0>(mincut);
            if(abs(round(weight)-weight) < pow(10,-10)) {
                weight = round(weight);
            }
        }
        if(!is_binary(x)) {
        // Branch
            auto frac_x = next(x.begin(),frac);
            GRBLinExpr branch_expr = frac_x->second;
            pair<GRBLinExpr,GRBLinExpr> branch0 = pair<GRBLinExpr,GRBLinExpr>(branch_expr,0);
            pair<GRBLinExpr,GRBLinExpr> branch1 = pair<GRBLinExpr,GRBLinExpr>(branch_expr,1);
            tsp(model, x, branch1, best_x, best_cost, n, m, nodes);
            tsp(model, x, branch0, best_x, best_cost, n, m, nodes);
            model.remove(constraint);
        } else {
        // New lower bound
            best_cost = C;
            best_x = x_val;
            model.remove(constraint);
        }
    }
    return best_cost;
}

double tsp(GRBModel& model, map<pair<int,int>,GRBVar>& x, pair<GRBLinExpr,GRBLinExpr>& constr, 
           vector<double>& best_x, double& best_cost, int n, int m, vector<int>& nodes) {
    auto constraint = model.addConstr(constr.first == constr.second);
    model.optimize();
    vector<double> x_val = vector<double>(m,0);
    map<pair<int,int>,double> lp_edge_map;
    
    if(model.get(GRB_IntAttr_Status) == 3) {
        model.remove(constraint);
        return best_cost;
    }

    double C = model.get(GRB_DoubleAttr_ObjVal);
    if(C >= best_cost) {
        model.remove(constraint);
        return best_cost;
    } else {
        for(int i=0; i<m; i++) {
            x_val[i] = (next(x.begin(),i)->second).get(GRB_DoubleAttr_X);
            lp_edge_map[next(x.begin(),i)->first] = x_val[i];
        }

        unsigned int frac = -1;
        double frac_diff = 1.0;
        for(int i=0; i<m; i++) {
            if(abs(x_val[i] - 0.5) < frac_diff) {
                frac_diff = abs(x_val[i] - 0.5);
                frac = i;
            }
        }
        auto mincut = min_cut(nodes,lp_edge_map,rand()%n);
        double weight = get<0>(mincut);
        if(abs(round(weight)-weight) < pow(10,-10)) {
            weight = round(weight);
        }

        while( weight < 2 ) {
        // Subtour
            GRBLinExpr st_expr = 0;
            auto mincut_edges = get<1>(mincut);
            mincut_edges = get<1>(mincut);
            for(int i=0; i<mincut_edges.size(); i++) {
                st_expr += x[mincut_edges[i]];
            }
            model.addConstr(st_expr >= 2);
            model.optimize();
            if(model.get(GRB_IntAttr_Status) == 3) {
                model.remove(constraint);
                return best_cost;
            }
            C = model.get(GRB_DoubleAttr_ObjVal);
            if(C >= best_cost) {
                model.remove(constraint);
                return best_cost;
            }
            for(int i=0; i<m; i++) {
                x_val[i] = (next(x.begin(),i)->second).get(GRB_DoubleAttr_X);
                lp_edge_map[next(x.begin(),i)->first] = x_val[i];
            }
            frac = -1;
            frac_diff = 1.0;
            for(int i=0; i<m; i++) {
                if( (abs(round(x_val[i]) - x_val[i]) > pow(10,-10)) ) {
                    if(abs(x_val[i] - 0.5) < frac_diff) {
                        frac_diff = abs(x_val[i] - 0.5);
                        frac = i;
                    }
                }
            }
            mincut = min_cut(nodes,lp_edge_map,rand()%n);
            weight = get<0>(mincut);
            if(abs(round(weight)-weight) < pow(10,-10)) {
                weight = round(weight);
            }
        }
        if(!is_binary(x)) {
        // Branch
            auto frac_x = next(x.begin(),frac);
            GRBLinExpr branch_expr = frac_x->second;
            pair<GRBLinExpr,GRBLinExpr> branch1 = pair<GRBLinExpr,GRBLinExpr>(branch_expr,1);
            pair<GRBLinExpr,GRBLinExpr> branch0 = pair<GRBLinExpr,GRBLinExpr>(branch_expr,0);
            tsp(model, x, branch1, best_x, best_cost, n, m, nodes);
            tsp(model, x, branch0, best_x, best_cost, n, m, nodes);
            model.remove(constraint);
        } else {
        // New lower bound
            best_cost = C;
            best_x = x_val;
            model.remove(constraint);
        }
    }
    return best_cost;
}

int main(int argc, char* argv[]) {
    ifstream f_in;
    ofstream f_out;
    int n; // # of nodes
    int m; // # of edges
    double best_obj = numeric_limits<double>::max();
    map<pair<int,int>,double> edge_map;
    map<pair<int,int>,double> lp_edge_map;
    vector<int> nodes;
    tuple<double,vector<pair<int,int> > > mincut;
 
    if(argc < 2) {
        cerr << "Usage: tsp input_file.txt" << endl;
        return 0;
    }
    f_in.open(argv[1]);
    if(!f_in.is_open()) {
        cerr << "Could not open tsp file!" << endl;
        return -1;
    }
    srand(time(NULL));
// Read node/edge sizes
    f_in >> n >> m;
// Read edges
    for(int i=0; i<m; i++) {
        int start, end, w;
        f_in >> start >> end >> w;
        edge_map[pair<int,int>(start,end)] = w;
        lp_edge_map[pair<int,int>(start,end)] = 1;
    }
// Set Nodes
    for(int i=0; i<n; i++) {
        nodes.push_back(i);
    }

    try {
        // Declare the GRB environment
        GRBEnv env = GRBEnv();
        // set the model within the GRB environment
        GRBModel model= GRBModel(env);
        // Set output to log
        model.set(GRB_IntParam_LogToConsole, 0);
        // Variables x_e
        map<pair<int,int>,GRBVar> x;
        vector<double> x_val = vector<double>(m,0);
        vector<double> best_x = vector<double>(m,0);
        for(auto it=edge_map.begin(); it!=edge_map.end(); ++it) {
            x[it->first] = model.addVar(0,1,0,GRB_CONTINUOUS);
        }
        model.update();
        // Objective
        GRBLinExpr master_obj = 0;
        for(auto it=x.begin(); it!=x.end(); ++it) {
            master_obj += (it->second)*edge_map[it->first];
        }
        model.setObjective(master_obj, GRB_MINIMIZE);
        // Degree Constraints
        vector<GRBLinExpr> vert_constr = vector<GRBLinExpr>(n,0);
        vector<GRBConstr> vertex_constraints = vector<GRBConstr>(n);
        for(auto it=x.begin(); it!=x.end(); ++it) {
            vert_constr[(it->first).first] += it->second;
            vert_constr[(it->first).second] += it->second;
        }
        for(int i=0; i<n; i++) {
            vertex_constraints[i] = model.addConstr(vert_constr[i] == 2);
        }
        // LP Stack
        stack<pair<GRBLinExpr,GRBLinExpr> > S;
        S.push(pair<GRBLinExpr,GRBLinExpr>(2,2));

        stack<GRBConstr> S_constr;

        int iterations = 0;
        int subtour_iter = 0;
        double greedy;
        for(int i=0; i<n; i++) {
            greedy = nearest_neighbor(nodes, edge_map, i);
            if(greedy < best_obj)
                best_obj = greedy;
        }
        cout << "Greedy upperbound Cost: " << best_obj << endl;

#ifdef STACK
        auto begin = chrono::high_resolution_clock::now();
        while(S.size() != 0) {
            model.update();

            S_constr.push(model.addConstr(S.top().first == S.top().second));

            model.optimize();

            if(model.get(GRB_IntAttr_Status) == 3) {
                iterations++;
                while( (S.top().second.getValue() == 0) ) {
                    model.remove(S_constr.top());
                    S.pop();
                    S_constr.pop();
                }
                model.remove(S_constr.top());
                S.pop();
                S_constr.pop();

            } else {
                double C = model.get(GRB_DoubleAttr_ObjVal);
                if(C >= best_obj) {
            iterations++;
                    // Remove last constraint
                    while( (S.top().second.getValue() == 0) ) {
                        model.remove(S_constr.top());
                        S.pop();
                        S_constr.pop();
                    }
                    model.remove(S_constr.top());
                    S.pop();
                    S_constr.pop();
                } else {
                    for(int i=0; i<m; i++) {
                        x_val[i] = (next(x.begin(),i)->second).get(GRB_DoubleAttr_X);
                        lp_edge_map[next(x.begin(),i)->first] = x_val[i];
                    }

                    unsigned int frac = -1;
                    double frac_diff = 1.0;
                    for(int i=0; i<m; i++) {
                        if(abs(x_val[i] - 0.5) < frac_diff) {
                            frac_diff = abs(x_val[i] - 0.5);
                            frac = i;
                        }
                    }
                    mincut = min_cut(nodes,lp_edge_map,next(x.begin(),frac)->first.first);

                // Subtour elimination constraints
                    double weight = get<0>(mincut);
                    if(abs(round(weight)-weight) < pow(10,-7)) {
                        weight = round(weight);
                    }
                    if( weight < 2) {
                        GRBLinExpr st_expr = 0;
                        auto mincut_edges = get<1>(mincut);
                        for(int i=0; i<mincut_edges.size(); i++) {
                            st_expr += x[mincut_edges[i]];
                        }
                        model.addConstr(st_expr >= 2);
                        model.remove(S_constr.top());
                        model.write("current_lp.lp");
                        subtour_iter++;
                        S_constr.pop();
                    } else if(!is_binary(x)) {
            iterations++;
                    // Find fractional x
                        auto frac_x = next(x.begin(),frac);
                        GRBLinExpr branch_expr = frac_x->second;
                    // Remove last constraint from model
                        S.push(pair<GRBLinExpr,GRBLinExpr>(branch_expr,0));
                        S.push(pair<GRBLinExpr,GRBLinExpr>(branch_expr,1));
                    } else {
            iterations++;
                        best_obj = C;
                        best_x = x_val;
                        while( (S.top().second.getValue() == 0)  ) {
                            model.remove(S_constr.top());
                            S.pop();
                            S_constr.pop();
                        }
                        model.remove(S_constr.top());
                        S.pop();
                        S_constr.pop();
                    }               
                }
            }
        }
        auto end = chrono::high_resolution_clock::now();
#endif
cout << "Best Cost from greedy algorithm: " << best_obj << endl;
#ifdef RECURSIVE
        auto begin = chrono::high_resolution_clock::now();
	double test_obj = tsp(model, x, S.top(), best_x, best_obj, n, m, nodes);
        auto end = chrono::high_resolution_clock::now();
#endif
        cout << "High Resolution Time: " << chrono::duration_cast<chrono::nanoseconds>(end-begin).count() << "ns" << endl;
        cout << "High Resolution Time (s): " << setprecision(12) << (double)(chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/pow(10,9)) << "s" << endl;
        cout << "Cost of best tour: " << best_obj << endl;       
#ifdef STACK
        cout << "Iterations: " << iterations << endl;
        cout << "Subtour Iterations: " << subtour_iter << endl;
#endif

// Output Results
        if(argc == 3) {
            f_out.open(argv[2]);
        } else {
            f_out.open("test_output.out");
        }
        auto edge_map_iter = edge_map.begin();
        for(int i=0; i<m; i++) {
            if(best_x[i] == 1) {
            	f_out << (edge_map_iter->first).first << " " << (edge_map_iter->first).second << " " << best_x[i] << " //for edge(" << i << ") of tour" << endl;
            }
            ++edge_map_iter;
        }
        f_out << "The cost of the best tour is: "<< best_obj << endl;
        f_out << "High Resolution Time: " << chrono::duration_cast<chrono::nanoseconds>(end-begin).count() << "ns" << endl;
        f_out << "High Resolution Time (s): " << setprecision(12) << (double)(chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/pow(10,9)) << "s" << endl;
        f_out << "Iterations: " << iterations << endl;
        f_out << "Subtour Iterations: " << subtour_iter << endl;
        f_out.close();

    } catch(GRBException e) {
        cout << "Gurobi Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during optimization" << endl;
    }
    return 1;
}

