/* PRL 2018 - Projekt 3 - Preorder prechod stromom
 * Autor: Marian Orszagh (xorsza00)
 *
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <cmath>

using namespace std;

#define RANK rank+1
// Ak je zadefinovana, normalny vystup je potlaceny a na vystup je vypisany
// cas behu
// #define TIME true


/* Vrchol zisti, ci ma laveho potomka.
 */
bool i_have_left_child(int node, int input_len) {
    return ((node+1)*2 <= input_len);
}


/* Vrchol zisti, ci ma praveho potomka.
 */
bool i_have_right_child(int node, int input_len) {
    return ((node+1)*2+1 <= input_len);
}


/* Vrchol zisti, ci ma predchodcu.
 */
bool i_have_parent(int node) {
    return (node > 0);
}


pair<int, int> rev_pair(pair<int, int> p) {
    return make_pair(p.second, p.first);
}


vector<pair<int, int>> tree_to_graph(int rank, int length) {
    vector<pair<int, int>> output_graph;
    if (i_have_left_child(rank, length)) {
        output_graph.push_back(make_pair(rank, rank*2));
        output_graph.push_back(make_pair(rank*2, rank));
    }
    if (i_have_right_child(rank, length)) {
        output_graph.push_back(make_pair(rank, rank*2+1));
        output_graph.push_back(make_pair(rank*2+1, rank));
    }
    return output_graph;
}


/* Pre hranu sa vypocita adjacency list, v ktorom sa nachadza.
 */
vector<int> calculate_adj_list(pair<int, int> edge, int length) {
    vector<int> adj_list;

    // Vyskyt laveho potomka
    if (i_have_left_child(edge.second, length)) {
        adj_list.push_back(edge.second * 4);
        adj_list.push_back(edge.second * 4 + 1);
    }
    // Vyskyt praveho potomka
    if (i_have_right_child(edge.second, length)) {
        adj_list.push_back(edge.second * 4 + 2);
        adj_list.push_back(edge.second * 4 + 3);
    }
    // Vyskyt rodica (v podstate vsetko okrem poslednej hrany etour)
    if (i_have_parent(edge.second)) {
        adj_list.push_back((edge.second-1) * 2 + 1);
        adj_list.push_back((edge.second-1) * 2);
    }

    return adj_list;
}


/* Vypocet nasledujucej hrany z korespondujuceho adjacency listu.
 */
int get_succ_from_adjlist(vector<int> adj_list, int rev_pid) {
    int my_succ;
    for (int i=0; i < adj_list.size(); i=i+2) {
        if (adj_list[i] == rev_pid) {
            if(i+2 < adj_list.size()){
                my_succ = adj_list[i+2];
            }
            else my_succ = adj_list[0];
        }
    }
    return my_succ;
}


int main(int argc, char *argv[])
{
    // Init
    int proc_count;
    int pid, rev_pid;
    int my_succ, my_pred;
    MPI_Status stat;

    // Meranie casu
    #ifdef TIME
        double start, finish;
    #endif

    // Vstupny retazec
    string input = argv[1];
    int length = input.size();

    if (length == 1) {
        cout << input[0] << endl;
        return 0;
    }

    // MPI init
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    #ifdef TIME
        if (!pid) start=MPI_Wtime();
    #endif

    pair<int, int> my_edge;
    pair<int, int> rev_edge;

    // Kazdy procesor spocita hranu svoju, reverznu k svojej a pid procesoru
    // reverznej hrany.
    if (pid%2) {
        my_edge = make_pair(pid/2+1, pid/4);
        rev_edge = rev_pair(my_edge);
        rev_pid = pid - 1;
    }
    else {
        my_edge = make_pair(pid/4, pid/2+1);
        rev_edge = rev_pair(my_edge);
        rev_pid = pid + 1;
    }

    // Pre hranu sa spocita adjacency list, v ktorom sa nachadza
    vector<int> adj_list = calculate_adj_list(my_edge, length);

    // Hrana zisti, kto je jej naslednikom
    my_succ = get_succ_from_adjlist(adj_list, rev_pid);
    // cout<<pid<<"->"<<my_succ<<endl;

    // Vektor naslednikov - cyklicky
    vector<int> succ(proc_count);
    succ[pid] = my_succ;
    MPI_Allgather(&succ[pid], 1, MPI_INT, &succ[0], 1, MPI_INT,
                  MPI_COMM_WORLD);


    // Svojmu naslednikovi poslem svoje id a prijmem id svojho predchodcu
    MPI_Send(&pid, 1, MPI_INT, my_succ, 0, MPI_COMM_WORLD);
    MPI_Recv(&my_pred, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);

    MPI_Barrier(MPI_COMM_WORLD);

    // Suffix sum
    // Svojemu succ poslem svojho predchodzu a moje id (DEBUG).
    // (pred, my_id)
    pair<int, int> transport_pair;

    int my_weight;
    // Posledna hrana musi byt spatna, nastavime jej -1, aby sa s dal lepsie
    // kontrolovat koniec zoznamu
    if (my_succ == 0) {
        my_weight = 0;
        my_succ = -1;
    }
    else my_weight = (my_edge.first < my_edge.second);
    // Prva hrana nema predchodcu a tak nebude ocakavat ziadnu spravu, mohlo
    // by dojst k zaseknutiu.
    if (pid == 0) my_pred = -1;

    // cout << pid+1 << "->" << my_weight<<" ";

    for (int n=0; n <= ceil(log2(length)); n++) {
        // cout<<pid+1<<endl;
        if (my_succ != -1) {
            // Poslem svojemu succ id svojeho predchodcu v pripade, ze nie som
            // na konci.
            transport_pair.first = my_pred;
            transport_pair.second = pid;

            MPI_Send(&transport_pair, 2, MPI_INT, my_succ, 0, MPI_COMM_WORLD);
            // cout<<pid<<" requested from  "<<my_succ<<" | " <<my_pred<<endl;
        }
        if (my_pred != -1) {
            // V pripade, ze mam co prijimat, prijimam request.
            MPI_Recv(&transport_pair, 2, MPI_INT, my_pred, 0, MPI_COMM_WORLD,
                     &stat);

            // Zaznamenam si, ktory pc mi poslal request, aby som mu mohol
            // vratit odpoved.
            int pred_id = transport_pair.second;
            // Zmenim id svojho predchodcu.
            my_pred = transport_pair.first;

            // cout<<pid+1<<" sent to "<<pred_id+1<<endl;

            // Prvou polozkou bude moj succ.
            transport_pair.first = my_succ;
            // Druhou moja hodnota.
            transport_pair.second = my_weight;
            // Poslem naspat.
            MPI_Send(&transport_pair, 2, MPI_INT, pred_id, 0, MPI_COMM_WORLD);
        }
        if (my_succ != -1) {
            MPI_Recv(&transport_pair, 2, MPI_INT, my_succ, 0,
                     MPI_COMM_WORLD, &stat);
            // cout<<pid+1 << " received from " << my_succ+1<<endl;

            // cout<<pid+1 << ": " << my_succ+1<<endl;
            // Zmenim si naslednika a urobim vypocet weight.
            my_succ = transport_pair.first;
            my_weight = my_weight + transport_pair.second;
        }
        // cout<<pid+1 << " got through" <<endl;

        MPI_Barrier(MPI_COMM_WORLD);
        // cout<<"----------------"<<endl;
    }
    // cout<<pid+1<<":"<<my_weight<<endl;

    // Korekcia
    // Dvojica, kde prvy prvok je poradie a druhy je vrchol
    pair<int, int> preorder(-2,-2);
    vector<pair<int, int>> preorder_whole(proc_count);
    if (my_edge.first < my_edge.second) {
        preorder.second = my_edge.second;
        preorder.first = proc_count - my_weight + 1;
    }
    // 1ka je urcite spatna hrana
    if (pid == 1) {
        preorder.first = -1;
        preorder.second = 0;
    }
    preorder_whole[pid] = preorder;
    MPI_Allgather(&preorder_whole[pid], 2, MPI_INT, &preorder_whole[0], 2,
                  MPI_INT, MPI_COMM_WORLD);

    // Distribucia vysledku
    preorder_whole[pid] = preorder;
    if (pid == 0) {
        // for (int i=0; i<preorder_whole.size(); i++) {
        //     cout<<preorder_whole[i].first<<":"<<preorder_whole[i].second+1<<endl;
        // }
        #ifndef TIME
        sort(preorder_whole.begin(), preorder_whole.end());
        for (int i=preorder_whole.size()-length;
             i<preorder_whole.size();
             i++) {
            cout<<input[preorder_whole[i].second];
        }
        cout<<endl;
        #else
        finish=MPI_Wtime();
        cout<<std::setprecision(40)<<(finish-start)<<endl;
        #endif
    }

    MPI_Finalize();
    return 0;
}
