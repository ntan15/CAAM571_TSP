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
#include <stack>
#include <functional>
#include "gurobi_c++.h"

using namespace std;

unsigned int hash_test(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

int contains_val(vector< vector<int> >& v, int val);
tuple<int,vector<pair<int,int> > > min_cut_phase(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e, int a);
tuple<int,vector<pair<int,int> > > min_cut(vector<int>& v, map<pair<int,int>,int>& e, int a);
int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e);
bool is_binary(map<pair<int,int>,GRBVar>& x, int m);


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
                            weight += iter->second;
                            edges.push_back(make_pair(B[k][l],A[i][j]));
                        }
                    }
                }
            }
        }
    }
    return make_tuple(weight,edges);
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

int most_tight(vector< vector<int> >& A, vector< vector<int> >& B, map<pair<int,int>,int>& e)
{
    int idx = -1;
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

bool is_binary(map<pair<int,int>,GRBVar>& x) {
    for(auto it=x.begin(); it!=x.end(); ++it) {
        double xval = (it->second).get(GRB_DoubleAttr_X);
        if(floor(xval) != xval) {
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    ifstream f_in;
    ofstream f_out;
    int n; // # of nodes
    int m; // # of edges
    double best_obj = numeric_limits<double>::max();
    map<pair<int,int>,int> edge_map, lp_edge_map;
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
for(auto it=edge_map.begin(); it!=edge_map.end(); ++it) {
cout << (it->first).first << " " << (it->first).second << " " << it->second << endl;
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
    	model.write("initial.lp");
        // Last added branch constraint
        GRBConstr constr;
        pair<GRBLinExpr,GRBLinExpr> constr_expr;
        // LP Stack
        stack<pair<GRBLinExpr,GRBLinExpr> > S;
        GRBLinExpr initial_constr = 0;
        S.push(pair<GRBLinExpr,GRBLinExpr>(initial_constr,initial_constr));
clock_t timer = clock();
int iterations = 0;
int subtour_iter = 0;
        while(S.size() != 0) {
iterations++;
model.update();
//cout << "Getting next branch constraint" << endl;
            constr_expr = S.top();
//cout << "Adding branch constraint to model" << endl;
            constr = model.addConstr(constr_expr.first == constr_expr.second);
            model.optimize();
//cout << "Model Optimized, Retrieving Objective Value" << endl;
            double C = model.get(GRB_DoubleAttr_ObjVal);
            if(C >= best_obj) {
//cout << "Objective value is worse than current best" << endl;
//cout << "    Removing branch constraint from stack and moving to next constraint" << endl;
                // Remove last constraint
                model.remove(constr);
                S.pop();
                // Add next constraint
//                constr_expr = S.top();
//                constr = model.addConstr(constr_expr.first == constr_expr.second);
            } else {
//cout << "Objective value is better than current best" << endl;
                // Modify edge weights to correspond to solution of LP
//cout << "Retrieving solution" << endl;
/*
for(int i=0; i<m; i++) {
cout << " " << (next(x.begin(),i)->second).get(GRB_DoubleAttr_X);
}
cout << endl;
*/
auto testx = model.getVars();
for(int i=0; i<m; i++) {
cout << testx[i].get(GRB_DoubleAttr_X) << " ";
}
cout << endl;

                for(int i=0; i<m; i++) {
                    x_val[i] = (next(x.begin(),i)->second).get(GRB_DoubleAttr_X);
cout << " " << x_val[i];
                    lp_edge_map[next(x.begin(),i)->first] = floor(x_val[i]);
                }
cout << endl;
/*
for(int i=0; i<m; i++) {
cout << " " << x_val[i];
}
cout << endl;
*/
//cout << "Generating Mincut" << endl;                
                // Generate mincut based on LP solution
                int vert = rand() % n;
cout << "Starting vertex for mincut: " << vert << endl;
                tuple<int,vector<pair<int,int> > > mincut = min_cut(nodes,lp_edge_map,vert);
//                tuple<int,vector<pair<int,int> > > mincut = min_cut(nodes,edge_map,vert);

vector<pair<int,int> > min_cut_edges = get<1>(mincut);
for(int i=0; i<min_cut_edges.size(); i++) {
    cout << min_cut_edges[i].first << " " << min_cut_edges[i].second << endl;
}

cout << "Mincut weight: " << get<0>(mincut) << endl;                
                // Subtour elimination constraints
                if(get<0>(mincut) < 2) {
//cout << "Mincut violated, adding subtour constraint" << endl;
                    GRBLinExpr st_expr = 0;
                    auto mincut_edges = get<1>(mincut);
                    for(int i=0; i<mincut_edges.size(); i++) {
                        st_expr += x[mincut_edges[i]];
                    }
                    model.addConstr(st_expr >= 2);
                    model.remove(constr);
subtour_iter++;
                } else if(!is_binary(x)) {
//cout << "Solution satisfies subtour but is fractional" << endl;
                    // Find fractional x
                    unsigned int frac = -1;
                    for(frac=0; frac<m; frac++) {
                        if(floor(x_val[frac]) != x_val[frac]) {
                            break;
                        }
                    }
                    auto frac_x = next(x.begin(),frac);
                    GRBLinExpr branch_expr = frac_x->second;
                    // Remove last constraint from model
//cout << "    Removing last branch constraint from model" << endl;                    
                    model.remove(constr);
//cout << "    Removing last branch constraint from stack" << endl;                    
                    S.pop();
//cout << "    Adding 2 branch constraints" << endl;
                    S.push(pair<GRBLinExpr,GRBLinExpr>(branch_expr,0));
                    S.push(pair<GRBLinExpr,GRBLinExpr>(branch_expr,1));
                    // Add next constraint to model
//                    constr_expr = S.top();
//                    constr = model.addConstr(constr_expr.first == constr_expr.second);
                } else {
//cout << "Solution satisfies subtour, is fractional, and is better than current best" << endl;
                    best_obj = C;
                    best_x = x_val;
/*
for(int i=0; i<m; i++) {
cout << " " << best_x[i];
}
cout << endl;
*/
//cout << "    Removing last branch constraint from model" << endl;                    
                    // Remove last constraint from model
                    model.remove(constr);
//cout << "    Removing last branch constraint from stack" << endl;                    
                    S.pop();
                    // Add next constraint to model
//cout << "   Getting next branch constraint" << endl;
//                    constr_expr = S.top();
//cout << "   Adding next branch constraint from model" << endl;
//                    constr = model.addConstr(constr_expr.first == constr_expr.second);
                }               
            }
//        } while(S.size() != 0);
        }
timer = clock() - timer;
        cout << "Cost of best tour: " << best_obj << endl;       
cout.precision(numeric_limits<double>::max_digits10);
cout << "Algorithm Time: " << scientific << ((double)timer)/(double)CLOCKS_PER_SEC << endl;
cout << "Algorithm Clock Tics: " << timer << endl;
cout << "Iterations: " << iterations << endl;
cout << "Subtour Iterations: " << subtour_iter << endl;

// Output Results
        f_out.open("test_output.out");
        auto edge_map_iter = edge_map.begin();
        for(int i=0; i<m; i++) {
            if(best_x[i] == 1) {
            	f_out << (edge_map_iter->first).first << " " << (edge_map_iter->first).second << " " << best_x[i] << " //for edge(" << i << ") of tour" << endl;
            }
            ++edge_map_iter;
        }
        f_out << "The cost of the best tour is: "<< best_obj << "(the cost the best tour)" << endl;
        f_out.close();

    } catch(GRBException e) {
        cout << "Gurobi Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during optimization" << endl;
    }


clock_t timer = clock();
tuple<int,vector<pair<int,int> > > test_map = min_cut(nodes,edge_map,rand()%n);
timer = clock() - timer;
cout.precision(numeric_limits<double>::max_digits10);
cout << "Min Cut Algorithm Map Time: " << scientific << ((double)timer)/(double)CLOCKS_PER_SEC << endl;
cout << "Min Cut Algorithm Map Clock Tics: " << timer << endl;
int weight = get<0>(test_map);
vector<pair<int,int> > min_edges = get<1>(test_map);
cout << "Cost of Min cut: " << weight << endl;

cout << "Starting at vertex 0" << endl;
test_map = min_cut(nodes,edge_map,rand()%n);
weight = get<0>(test_map);
min_edges = get<1>(test_map);
cout << "Cost of Min cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

cout << "Starting at vertex 1" << endl;
test_map = min_cut(nodes,edge_map,1);
weight = get<0>(test_map);
min_edges = get<1>(test_map);
cout << "Cost of Min Cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

cout << "Starting at vertex 2" << endl;
test_map = min_cut(nodes,edge_map,2);
weight = get<0>(test_map);
min_edges = get<1>(test_map);
cout << "Cost of Min Cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

cout << "Starting at vertex 3" << endl;
test_map = min_cut(nodes,edge_map,3);
weight = get<0>(test_map);
min_edges = get<1>(test_map);
cout << "Cost of Min Cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

cout << "Starting at vertex 4" << endl;
test_map = min_cut(nodes,edge_map,4);
weight = get<0>(test_map);
min_edges = get<1>(test_map);
cout << "Cost of Min Cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

cout << "Starting at vertex 5" << endl;
test_map = min_cut(nodes,edge_map,5);
weight = get<0>(test_map);
min_edges = get<1>(test_map);
cout << "Cost of Min Cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

cout << "Starting at vertex 6" << endl;
test_map = min_cut(nodes,edge_map,6);
weight = get<0>(test_map);
min_edges = get<1>(test_map);
cout << "Cost of Min Cut: " << weight << endl;
for(int i=0; i<min_edges.size(); i++) {
    cout << get<0>(min_edges[i]) << " " << get<1>(min_edges[i]) << endl;
}

    return 1;
}

