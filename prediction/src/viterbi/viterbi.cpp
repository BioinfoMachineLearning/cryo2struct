#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <tuple> 
#include <set>


// g++ -fPIC -shared -o viterbi.so viterbi.cpp -O3
// // conda install -c conda-forge gxx  -> if older version of g++ is present

bool success = false;

std::vector<int> exclude_states_in_c; //global exclude states

// Function prototype for viterbi_main with C linkage specification
extern "C" void run_viterbi(std::vector<int> observation, int num_observations, std::vector<int> states, int num_states, std::vector<std::vector<double>> transition_matrix, std::vector<std::vector<double>> emission_matrix, std::vector<double> initial_matrix, std::vector<int> states_to_work_python);

// Implementation of my_cpp_function with C linkage specification
extern "C" void viterbi(std::vector<int> observation, int num_observations, std::vector<int> states, int num_states, std::vector<std::vector<double>> transition_matrix, std::vector<std::vector<double>> emission_matrix, std::vector<double> initial_matrix, std::vector<int> states_to_work_python) 
{
    // Run viterbi program
    std::vector<std::vector<double>> trellis(num_observations, std::vector<double>(num_states, 0));
    std::unordered_map<int, std::vector<int>> path;  //an empty unordered_map
    std::vector<int> states_to_work_with;
    for (int x=0; x<states_to_work_python.size(); x++)
    {
        if (std::find(exclude_states_in_c.begin(), exclude_states_in_c.end(), states_to_work_python[x]) == exclude_states_in_c.end())
        {
            states_to_work_with.push_back(states_to_work_python[x]);
        }
    }
    for (int state: states_to_work_with)
    {  
        trellis[0][state] = initial_matrix[state] + emission_matrix[state][observation[0]];
        path[state].push_back(state);
    }
    for (int observation_index=1; observation_index<num_observations; observation_index++)
    {
        std::unordered_map<int, std::vector<int>> new_path;
        for (int state: states_to_work_with)
        {
            double max_prob = -std::numeric_limits<float>::infinity();
            int possible_state = -1;
            for (int previous_state: states_to_work_with)
            {
                auto it = std::find(path[previous_state].begin(), path[previous_state].end(), state);
                if (it != path[previous_state].end())
                {
                    continue;
                }
                double prob = trellis[observation_index - 1][previous_state] + transition_matrix[previous_state][state] + emission_matrix[state][observation[observation_index]];
                if (prob > max_prob)
                {
                    max_prob = prob;
                    possible_state = previous_state;
                }
            }
            double probability = max_prob;
            if (possible_state == -1)
            {
                if (success == true){return;}
                std::cout<<" Successfully connected till " << observation_index << " of " << observation.size() << std::endl;
                std::cout<<" Running for fragments -- " <<std::endl;
                auto max_it = std::max_element(states_to_work_with.begin(), states_to_work_with.end(),
                               [&](const auto& state1, const auto& state2) {
                                   return trellis[observation_index - 1][state1] <
                                          trellis[observation_index - 1][state2];
                               });

                auto probability = trellis[observation_index - 1][*max_it];
                auto state = *max_it;
                exclude_states_in_c.insert(exclude_states_in_c.end(), path[state].begin(), path[state].end());
                std::vector<int> sub_observation(observation.begin() + observation_index, observation.end());
                run_viterbi(sub_observation, sub_observation.size(), states, num_states, transition_matrix, emission_matrix, initial_matrix, states_to_work_python);
            }
            trellis[observation_index][state] = probability;
            new_path[state] = path[possible_state];
            new_path[state].push_back(state);
        }
        path = new_path;
    }
    auto max_it = std::max_element(states_to_work_with.begin(), states_to_work_with.end(),
                    [&](const auto& state1, const auto& state2) {
                        return trellis[num_observations - 1][state1] <
                                trellis[num_observations - 1][state2];
                    });
    auto probability = trellis[num_observations - 1][*max_it];
    auto state = *max_it;   
    exclude_states_in_c.insert(exclude_states_in_c.end(), path[state].begin(), path[state].end());
    success = true;
    return;
    
}


extern "C" void run_viterbi(std::vector<int> observation, int num_observations, std::vector<int> states, int num_states, std::vector<std::vector<double>> transition_matrix, std::vector<std::vector<double>> emission_matrix, std::vector<double> initial_matrix, std::vector<int> states_to_work_python)
{
    viterbi(observation, num_observations, states, num_states, transition_matrix, emission_matrix, initial_matrix, states_to_work_python); 
    return;
    
}


// Implementation of viterbi_main with C linkage specification
extern "C" int* viterbi_main(int* obs, int num_observations, int num_states, double* transt, double* emiss, double* init, int* exclude_s, int exclude_s_len) 
{
    std::cout << "############################ Running Alignment Program (VITERBI) ############################" << std::endl;
    success = false;
    std::vector<int> observations(num_observations);
    std::vector<int> states(num_states);

    std::vector<std::vector<double>> transition_matrix(num_states, std::vector<double>(num_states));
    std::vector<std::vector<double>> emission_matrix(num_states, std::vector<double>(20));
    std::vector<double> initial_matrix(num_states);
    std::vector<int> exclude_stat(exclude_s_len);

    for (int i = 0; i<num_states; i++)
    {
        initial_matrix[i] = init[i];
        states[i] = i;
        for (int j = 0; j<num_states; j++)
        {
            transition_matrix[i][j] = transt[i*num_states + j];
        }
    }
    
    for (int i = 0; i<num_states; i++)
    {
        for (int j = 0; j<20; j++)
        {
            emission_matrix[i][j] = emiss[i*20 + j];
        }
    }

    for (int i = 0; i<exclude_s_len; i++)
    {
        exclude_stat[i] = exclude_s[i];
    }

    for (int i = 0; i<num_observations; i++)
    {
        observations[i] = obs[i];
    }

    std::vector<int> states_to_work_python;
    for (int x=0; x<num_states; x++)
    {
        if (std::find(exclude_stat.begin(), exclude_stat.end(), states[x]) == exclude_stat.end())
        {
            states_to_work_python.push_back(states[x]);
        }
    }
    
    exclude_states_in_c.clear();
    run_viterbi(observations, num_observations, states, num_states, transition_matrix, emission_matrix, initial_matrix, states_to_work_python);
 
    std::cout << "Length of exclude states " << exclude_states_in_c.size() << std::endl;
    std::set<int> s(exclude_states_in_c.begin(), exclude_states_in_c.end());
    std::cout << "Number of unique values: " << s.size() << std::endl;
    int * result = exclude_states_in_c.data();
    return result;
    
}