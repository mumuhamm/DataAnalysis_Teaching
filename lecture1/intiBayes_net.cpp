
/*
    

    The network contains 4 nodes and looks as follows:

        B     C
        \\   //
         \/ \/ 
           A
           ||
           \/
            D


    The probabilities of each node are summarized below.  (The probability
    of each node being 0 is not listed since it is just P(X=0) = 1-p(X=1) ) 

        p(B=1) = 0.01

        p(C=1) = 0.001

        p(A=1 | B=0, C=0) = 0.01  
        p(A=1 | B=0, C=1) = 0.5
        p(A=1 | B=1, C=0) = 0.9
        p(A=1 | B=1, C=1) = 0.99 

        p(D=1 | A=0) = 0.2 
        p(D=1 | A=1) = 0.5

*/


#include <dlib/bayes_utils.h>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <dlib/directed_graph.h>
#include <iostream>


using namespace dlib;
using namespace std;

// ----------------------------------------------------------------------------------------

int main()
{
    try
      {
        
        using namespace bayes_node_utils;
        directed_graph<bayes_node>::kernel_1a_c bn;
        enum nodes
        {
	  A = 0,
	  B = 1,
	  C = 2,
            D = 3
        };

      
        bn.set_number_of_nodes(4);
        bn.add_edge(A, D);
        bn.add_edge(B, A);
        bn.add_edge(C, A);


         
        set_node_num_values(bn, A, 2);
        set_node_num_values(bn, B, 2);
        set_node_num_values(bn, C, 2);
        set_node_num_values(bn, D, 2);

        assignment parent_state;
        
        set_node_probability(bn, B, 1, parent_state, 0.01);
       
        set_node_probability(bn, B, 0, parent_state, 1-0.01);


       
        set_node_probability(bn, C, 1, parent_state, 0.001);
        
        set_node_probability(bn, C, 0, parent_state, 1-0.001);


        parent_state.add(B, 1);
        parent_state.add(C, 1);
      
        set_node_probability(bn, A, 1, parent_state, 0.99);
        // Here we specify that p(A=0 | B=1, C=1) = 1-0.99 
        set_node_probability(bn, A, 0, parent_state, 1-0.99);

        // Here we use the [] notation because B and C have already
        // been added into parent state.  
        parent_state[B] = 1;
        parent_state[C] = 0;
        // Here we specify that p(A=1 | B=1, C=0) = 0.9 
        set_node_probability(bn, A, 1, parent_state, 0.9);
        set_node_probability(bn, A, 0, parent_state, 1-0.9);

        parent_state[B] = 0;
        parent_state[C] = 1;
        // Here we specify that p(A=1 | B=0, C=1) = 0.5 
        set_node_probability(bn, A, 1, parent_state, 0.5);
        set_node_probability(bn, A, 0, parent_state, 1-0.5);

        parent_state[B] = 0;
        parent_state[C] = 0;
        // Here we specify that p(A=1 | B=0, C=0) = 0.01 
        set_node_probability(bn, A, 1, parent_state, 0.01);
        set_node_probability(bn, A, 0, parent_state, 1-0.01);

        parent_state.clear();
        parent_state.add(A,1);
        // Here we specify that p(D=1 | A=1) = 0.5 
        set_node_probability(bn, D, 1, parent_state, 0.5);
        set_node_probability(bn, D, 0, parent_state, 1-0.5);

        parent_state[A] = 0;
        // Here we specify that p(D=1 | A=0) = 0.2 
        set_node_probability(bn, D, 1, parent_state, 0.2);
        set_node_probability(bn, D, 0, parent_state, 1-0.2);


        typedef dlib::set<unsigned long>::compare_1b_c set_type;
        typedef graph<set_type, set_type>::kernel_1a_c join_tree_type;
        join_tree_type join_tree;

     
        create_moral_graph(bn, join_tree);
        create_join_tree(join_tree, join_tree);

       
        bayesian_network_join_tree solution(bn, join_tree);


        // now print out the probabilities for each node
        cout << "Using the join tree algorithm:\n";
        cout << "p(A=1) = " << solution.probability(A)(1) << endl;
        cout << "p(A=0) = " << solution.probability(A)(0) << endl;
        cout << "p(B=1) = " << solution.probability(B)(1) << endl;
        cout << "p(B=0) = " << solution.probability(B)(0) << endl;
        cout << "p(C=1) = " << solution.probability(C)(1) << endl;
        cout << "p(C=0) = " << solution.probability(C)(0) << endl;
        cout << "p(D=1) = " << solution.probability(D)(1) << endl;
        cout << "p(D=0) = " << solution.probability(D)(0) << endl;
        cout << "\n\n\n";


        set_node_value(bn, C, 1);
        set_node_as_evidence(bn, C);

       
        bayesian_network_join_tree solution_with_evidence(bn, join_tree);

        // now print out the probabilities for each node
        cout << "Using the join tree algorithm:\n";
        cout << "p(A=1 | C=1) = " << solution_with_evidence.probability(A)(1) << endl;
        cout << "p(A=0 | C=1) = " << solution_with_evidence.probability(A)(0) << endl;
        cout << "p(B=1 | C=1) = " << solution_with_evidence.probability(B)(1) << endl;
        cout << "p(B=0 | C=1) = " << solution_with_evidence.probability(B)(0) << endl;
        cout << "p(C=1 | C=1) = " << solution_with_evidence.probability(C)(1) << endl;
        cout << "p(C=0 | C=1) = " << solution_with_evidence.probability(C)(0) << endl;
        cout << "p(D=1 | C=1) = " << solution_with_evidence.probability(D)(1) << endl;
        cout << "p(D=0 | C=1) = " << solution_with_evidence.probability(D)(0) << endl;
        cout << "\n\n\n";

       
        set_node_value(bn, A, 0);
        set_node_value(bn, B, 0);
        set_node_value(bn, D, 0);

        bayesian_network_gibbs_sampler sampler;


        unsigned long A_count = 0;
        unsigned long B_count = 0;
        unsigned long C_count = 0;
        unsigned long D_count = 0;

       
        const long rounds = 2000;
        for (long i = 0; i < rounds; ++i)
	  {
            sampler.sample_graph(bn);

            if (node_value(bn, A) == 1)
	      ++A_count;
            if (node_value(bn, B) == 1)
	      ++B_count;
            if (node_value(bn, C) == 1)
	      ++C_count;
            if (node_value(bn, D) == 1)
	      ++D_count;
	  }

        cout << "Using the approximate Gibbs Sampler algorithm:\n";
        cout << "p(A=1 | C=1) = " << (double)A_count/(double)rounds << endl;
        cout << "p(B=1 | C=1) = " << (double)B_count/(double)rounds << endl;
        cout << "p(C=1 | C=1) = " << (double)C_count/(double)rounds << endl;
        cout << "p(D=1 | C=1) = " << (double)D_count/(double)rounds << endl;
      }
    catch (std::exception& e)
      {
        cout << "exception thrown: " << endl;
        cout << e.what() << endl;
        cout << "hit enter to terminate" << endl;
        cin.get();
      }
}



