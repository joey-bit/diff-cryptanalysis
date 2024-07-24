#include <iostream>
#include <map>
#include <array>
#include <random>
#include <vector>
#include <algorithm>
#include <list>
#include <unordered_set>
#include <cstdlib> // For std::rand and std::srand
#include <ctime>   // For std::time
#include <cmath>


//Function Prototypes
//------------------------------------Fundamental SPN functions-------------------------------------------------//
std::array<int, 4> s_box_transform(const std::array<int, 4>&);
std::array<int, 16> permutation_transform(const std::array<int, 16>&);
std::array<int, 16> array_xor_sixteen(const std::array<int, 16>&, const std::array<int, 16>&); //for AddRoundKey
std::array<int, 16> run_spn(const std::array<int, 16>&); //run full SPN
void key_gen(void);
void subkey_five_attack(void);

//------------------------------------Probability and Value Generation Functions-------------------------------//
std::array<std::array<int, 16>, 16> build_characteristic(void);
std::vector<std::pair<int, int>> generate_input_pairs(std::array<int, 16>);
std::pair<std::array<int, 16>,float> expected_output_diff(const std::array<int, 16>&);
std::pair<std::vector<std::pair<int, int>>, float> isolate_right_pairs(const std::vector<std::pair<int, int>>&, const std::array<int, 16>&);
std::array<std::pair<std::pair<int, int>, float>, 5> compute_best_differences(void);

//------------------------------------Helper Functions---------------------------------------------------------//
std::array<std::array<int, 4>, 4> split_array(const std::array<int, 16>&);
std::array<int, 4> dec_to_binary_four(int);
std::array<int, 16> dec_to_binary_sixteen(int);
int binary_to_dec_four(const std::array<int, 4>&);
int binary_to_dec_sixteen(const std::array<int, 16>&);
std::array<int, 4> array_xor_four(const std::array<int, 4>&, const std::array<int, 4>&);
std::array<int, 16> concatenate_arrays(const std::array<std::array<int, 4>, 4>&);
int vector_to_hex(const std::vector<int>&);
std::vector<int> hex_to_vector(int);
std::vector<size_t> getRandomIndices(size_t, size_t);

//S-Box Definition
const std::map<int, int> s_box = {
    {0x0, 0xE},
    {0x1, 0x4},
    {0x2, 0xD},
    {0x3, 0x1},
    {0x4, 0x2},
    {0x5, 0xF},
    {0x6, 0xB},
    {0x7, 0x8},
    {0x8, 0x3},
    {0x9, 0xA},
    {0xA, 0x6},
    {0xB, 0xC},
    {0xC, 0x5},
    {0xD, 0x9},
    {0xE, 0x0},
    {0xF, 0x7}
}; 

const std::map<int, int> inverse_s_box = {
    {0xE, 0x0},
    {0x4, 0x1},
    {0xD, 0x2},
    {0x1, 0x3},
    {0x2, 0x4},
    {0xF, 0x5},
    {0xB, 0x6},
    {0x8, 0x7},
    {0x3, 0x8},
    {0xA, 0x9},
    {0x6, 0xA},
    {0xC, 0xB},
    {0x5, 0xC},
    {0x9, 0xD},
    {0x0, 0xE},
    {0x7, 0xF}
};

//Permutation Definition
const std::map<int, int> permutation = {
    {0, 0},
    {1, 4},
    {2, 8},
    {3, 12},
    {4, 1},
    {5, 5},
    {6 ,9},
    {7, 13},
    {8, 2},
    {9, 6},
    {10, 10},
    {11, 14},
    {12, 3},
    {13, 7},
    {14, 11},
    {15, 15},
};

//Keys for the SPN
std::array<std::array<int, 16>, 5> keys;

//Probability Characteristic for the S-Box
std::array<std::array<int, 16>, 16> prob_characteristic;

int main() {
    key_gen();
    prob_characteristic = build_characteristic();
    subkey_five_attack();
    
    return 0;
}

//random number generator from CPP reference
void key_gen(void) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, 1);
    for(auto& key : keys) {
        for(auto& bit : key) {
            bit = distrib(gen);
        }
    }
    std::cout << "Keys generated: \n";
    int i = 1;
    for (const auto& key : keys) {
        std::cout << "Subkey " << i++ << ": ";
        for (const auto& s : key) std::cout << s << ' ';
        std::cout << '\n';
    }
    std::cout << '\n';
}

//Transform a 4 bit binary input to a 4 bit binary output using s-box
std::array<int, 4> s_box_transform(const std::array<int, 4>& input) {
    std::array<int, 4> output;
    int decimal_value = input[0]*8 + input[1]*4 + input[2]*2 + input[3];
    int temp = s_box.at(decimal_value);
    //Adapted from geeksforgeeks.org
    for(auto i=3; i>=0; i--) {
        output[i] =  temp % 2;
        temp = temp / 2;
    }
    // std::cout << "Transformed value \n";
    // for (const auto& bit : output) std::cout << bit << " ";
    // std::cout << "\n";
    return output;
}

//Transform a 16 bit binary input to a binary output using permutation
std::array<int, 16> permutation_transform(const std::array<int, 16>& input) {
    std::array<int, 16> output;
    for(auto i=0; i<16; i++) {
        output[permutation.at(i)] = input[i];
    }
    // std::cout << "Transformed value \n";
    // for (const auto& bit : output) std::cout << bit << " ";
    // std::cout << "\n";
    return output;
}

//Concatenate 4 4-bit arrays to a single 16-bit array
std::array<int, 16> concatenate_arrays(const std::array<std::array<int, 4>, 4>& input) {
    std::array<int, 16> output;
    int i = 0;
    for(const auto& arr : input) {
        for(const auto& bit : arr) {
            output[i++] = bit;
        }
    }
    // std::cout << "Concatenated value \n";
    // for (const auto& bit : output) std::cout << bit << " ";
    // std::cout << "\n";
    return output;
}

//Disjoin a 16-bit array to 4 4-bit arrays
std::array<std::array<int, 4>, 4> split_array(const std::array<int, 16>& input) {
    std::array<std::array<int, 4>, 4> output;
    int j = 0;
    for(auto i=0; i<16; i++) {
        if(i%4 == 0 && i != 0) j++;
        output[j][i%4] = input[i];
    }
    // std::cout << "Split value \n";
    // for (const auto& arr : output) {
    //     for (const auto& bit : arr) std::cout << bit << " ";
    //     std::cout << "\n";
    // }
    return output;
}

//Convert a decimal number to a 4-bit binary array
std::array<int, 4> dec_to_binary_four(int decimal) {
    std::array<int, 4> binary;
    for(auto i=3; i>=0; i--) {
        binary[i] =  decimal % 2;
        decimal = decimal / 2;
    }
    // std::cout << "Binary value: ";
    // for (const auto& bit : binary) std::cout << bit << " ";
    return binary;
}

std::array<int, 16> dec_to_binary_sixteen(int decimal) {
    std::array<int, 16> binary;
    for(auto i=15; i>=0; i--) {
        binary[i] =  decimal % 2;
        decimal = decimal / 2;
    }
    // std::cout << "Binary value: ";
    // for (const auto& bit : binary) std::cout << bit << " ";
    return binary;
}

//Convert a 4-bit binary array to a decimal number
int binary_to_dec_four(const std::array<int, 4>& binary) {
    int decimal = 8*binary[0] + 4*binary[1] + 2*binary[2] + binary[3];
    return decimal;
}

//Convert a 16-bit binary array to a decimal number
int binary_to_dec_sixteen(const std::array<int, 16>& binary) {
    int decimal = 0;
    for(auto i=0; i<16; i++) {
        decimal += binary[i] * (1 << (15-i));
    }
    return decimal;
}

//XOR two 4-bit arrays
std::array<int, 4> array_xor_four(const std::array<int, 4>& arr1, const std::array<int, 4>& arr2) {
    std::array<int, 4> output;
    for(auto i=0; i<4; i++) {
        output[i] = arr1[i] ^ arr2[i];
    }
    return output;
}

//XOR two 16-bit arrays
std::array<int, 16> array_xor_sixteen(const std::array<int, 16>& arr1, const std::array<int, 16>& arr2) {
    std::array<int, 16> output;
    for(auto i=0; i<16; i++) {
        output[i] = arr1[i] ^ arr2[i];
    }
    return output;
}

//Build difference characteristics for a given s-box
std::array<std::array<int, 16>, 16> build_characteristic(void) {
    std::array<std::array<int, 16>, 16> characteristic = {0};
    for(auto i=0; i<16; i++) {
        for(auto j=0; j<16; j++) {
            auto input_diff = j^i; //input difference delta_x
            auto y_1 = s_box_transform(dec_to_binary_four(i)); //output y'
            auto y_2 = s_box_transform(dec_to_binary_four(j)); //output y''
            auto output_diff = binary_to_dec_four(array_xor_four(y_1, y_2)); //get delta_y
            characteristic[input_diff][output_diff]++;
        }
    }
    // std::cout << "Characteristics: \n";
    // for (const auto& row : characteristic) {
    //     for (const auto& cell : row) std::cout << cell << " ";
    //     std::cout << "\n";
    // }
    return characteristic;
}

//Run the SPN for a given plaintext
std::array<int, 16> run_spn(const std::array<int, 16>& plaintext) {
    std::array<int, 16> output;
    std::array<int, 16> input = plaintext; //initialize input to plaintext
    for (auto i=0; i<3; i++) {
        input = array_xor_sixteen(input, keys[i]); //add round key
        auto split = split_array(input); //split input into 4 4-bit arrays
        for (auto& array: split) array = s_box_transform(array); //apply s-box transform to each array
        input = concatenate_arrays(split); //concatenate the 4 4-bit arrays
        input = permutation_transform(input); //apply permutation transform
    }
    input = array_xor_sixteen(input, keys[3]); //add final round key
    auto split = split_array(input); //split input into 4 4-bit arrays
    for (auto& array: split) array = s_box_transform(array); //apply s-box transform to each array
    input = concatenate_arrays(split); //concatenate the 4 4-bit arrays
    output = array_xor_sixteen(input, keys[4]); //add final round key
    return output;
}

//Generate a dataset of plaintext pairs for which the pairs have a given input difference
std::vector<std::pair<int, int>> generate_input_pairs(std::array<int, 16> input_diff) {
    std::vector<std::pair<int, int>> input_pairs;
    const int input_diff_decimal = binary_to_dec_sixteen(input_diff);
    //std::cout << "Input difference: " << input_diff_decimal << "\n";
    for(auto i=0; i<65536; i++) {   //iterate through all possible plaintext pairs
        for(auto j=0; j<65536; j++){
            auto current_diff = i^j;
            if(current_diff != input_diff_decimal) continue;   //skip pairs that do not have the input difference
            input_pairs.push_back(std::make_pair(i, j));
        }
    }
    //std::cout << "Number of pairs generated: " << input_pairs.size() << "\n";
    return input_pairs;
}

//Compute the expected output difference for a given input difference, and it's probability
//Feed the input difference through the S box and permutation to get the output difference
std::pair<std::array<int, 16>,float> expected_output_diff(const std::array<int, 16>& input_diff) {
    auto input_arrays = split_array(input_diff); //split the input difference into 4 4-bit arrays
    float probability = 1; //initialize probability to 1

    //Compute output difference for each sub s-box
    for(auto i = 0; i<3; i++){
        for (auto& array : input_arrays) {
            if(binary_to_dec_four(array) == 0) continue; //skip arrays with 0 input difference
            //This section supported by GPT
            auto u_in = binary_to_dec_four(array); //get input difference
            auto max_prob_it = std::max_element(prob_characteristic[u_in].begin(),prob_characteristic[u_in].end()); //get the position of the maximum probability for the input difference
            probability *= float(*max_prob_it)/16; //multiply the probability of the output difference
            int v_out = std::distance(prob_characteristic[u_in].begin(), max_prob_it); //get the index of the maximum probability (which is the output difference)
            array = dec_to_binary_four(v_out); //convert the output difference to binary
        }
        auto output_diff = concatenate_arrays(input_arrays); //concatenate the 4 4-bit arrays
        output_diff = permutation_transform(output_diff); //apply permutation transform
        input_arrays = split_array(output_diff); //split the output difference into 4 4-bit arrays
    }
    //printf("Expected Output Probability: %3f\n", probability*100);
    return std::make_pair(concatenate_arrays(input_arrays), probability);

}

//Compute the best 100 input differences for the R-1 round differential characteristic
std::array<std::pair<std::pair<int, int>, float>, 5> compute_best_differences(void) {
    std::array<std::pair<std::pair<int, int>, float>, 5> best_differences;
    std::array<std::pair<std::pair<int, int>, float>, 65536> difference_probabiliies;
    for(auto i=0; i<65536; i++) {
        auto input_diff = dec_to_binary_sixteen(i);
        auto diff_prob = expected_output_diff(input_diff);
        std::pair<std::pair<int, int>, float> out;
        out.first.first = i;
        out.first.second = binary_to_dec_sixteen(diff_prob.first);
        out.second = diff_prob.second;
        difference_probabiliies[i] = out; //outputs stored as ((in, out), prob)
    }
    std::sort(difference_probabiliies.begin(), difference_probabiliies.end(), [](auto& left, auto& right) {
        return left.second > right.second;
    });
    for(auto i=1; i<6; i++) { //skip the first difference because the probability is 1
        best_differences[i-1] = difference_probabiliies[i];
    }
    std::cout << "Best Differences for Given S-Box: \n";
    for (const auto& diff : best_differences) {
        std::cout << "Input Difference: " << diff.first.first << " Output Difference: " << diff.first.second << " Probability: " << diff.second*100 << "%\n";
    }
    printf("\n");
    return best_differences;
}

//Issue an attack on the 5th subkey using the best 5 best input differences
void subkey_five_attack(void) {
    auto best_differences = compute_best_differences(); //Generate 5 input differences that have the highest probabilities of having the differential characteristic
    std::array<std::array<int, 16>, 4> target_partial_value_counter; //Initialize the counter for partial subkeys
    for(auto& counts : target_partial_value_counter) counts.fill(0);
    std::vector<std::pair<int, int>> ciphertext_pairs; //Generate a dataset of ciphertext pairs for the given input difference TODO
    std::vector<int> non_zero_indices; //store the indices of the 4 arrays so we know which keys are being tested
    std::vector<int> u_one; //use these to get the difference at the inputs
    std::vector<int> u_two;
    std::vector<int> non_zero_differences;
    auto output_diff = split_array(dec_to_binary_sixteen(best_differences[0].first.second)); //Get the last round difference we're looking for
    std::cout<<"Ideal Input Difference in Hex: ";
    for (const auto& diff: split_array(dec_to_binary_sixteen(best_differences[0].first.first))) std::cout << binary_to_dec_four(diff) << " ";
    std::cout << '\n';
    std::cout<<"Ideal Output Difference in Hex: ";
    for (const auto& diff: output_diff) std::cout << binary_to_dec_four(diff) << " ";
    std::cout << '\n';
    auto input_pairs = generate_input_pairs(dec_to_binary_sixteen(best_differences[0].first.first)); //Generate the input pairs based on input difference
    printf("Number of plaintext pairs generated: %d\n", input_pairs.size());
    auto random_indices = getRandomIndices(5000, input_pairs.size()); //generate an array of 5000 indices (i.e. randomly select input pairs)
    printf("Number of random plaintext pairs being tested %d\n", random_indices.size());

    printf("Plaintext pairs, listing 10 of 5000\n");
    for(auto i=0; i<10; i++) printf("P1: %d, P2: %d\n", input_pairs[random_indices[i]].first, input_pairs[random_indices[i]].second);
    printf("\n");
    for (const auto& index: random_indices){ //Generate 5000 pairs of ciphertexts
        auto ciphertext_one = binary_to_dec_sixteen(run_spn(dec_to_binary_sixteen(input_pairs[index].first)));
        auto ciphertext_two = binary_to_dec_sixteen(run_spn(dec_to_binary_sixteen(input_pairs[index].second)));
        ciphertext_pairs.push_back(std::make_pair(ciphertext_one, ciphertext_two));
    }

    printf("Ciphertext pairs generated, listing 10 of 5000\n");
    for(auto i=0; i<10; i++) printf("C1: %d, C2: %d\n", ciphertext_pairs[i].first, ciphertext_pairs[i].second);
    printf("\n");

    for(auto i=0; i<4; i++){
        if(binary_to_dec_four(output_diff[i]) == 0) continue; //only test active S boxes for the final difference
        non_zero_indices.push_back(i);
        non_zero_differences.push_back(binary_to_dec_four(output_diff[i]));
    }

    //calculate the number of possible keys
    auto possible_keys = static_cast<int>(pow(16, non_zero_indices.size()));
    printf("Number of possible keys being tested per ciphertext pair: %d\n", possible_keys);

    auto count = 0;

    for (auto i=0; i<10; i++) {
        //turn each ciphertext into a 4x4bit array
        auto out_one = split_array(dec_to_binary_sixteen(ciphertext_pairs[i].first)); 
        auto out_two = split_array(dec_to_binary_sixteen(ciphertext_pairs[i].second));
        
        //isolate the arrays corresponding to non-zero differences
        for(const auto& index: non_zero_indices){
            u_one.push_back(binary_to_dec_four(out_one[index]));
            u_two.push_back(binary_to_dec_four(out_two[index]));
        }

        for(auto j=0; j<possible_keys; j++){
            auto in_one = hex_to_vector(vector_to_hex(u_one)^j);
            auto in_two = hex_to_vector(vector_to_hex(u_two)^j);
            for(auto& val : in_one) val = inverse_s_box.at(val);
            for(auto& val : in_two) val = inverse_s_box.at(val);

            auto in_diff = vector_to_hex(in_one)^vector_to_hex(in_two);
            auto right_pair = in_diff == vector_to_hex(non_zero_differences);
            if(right_pair){
                // auto current_indices = hex_to_vector(j);
                // printf("Key match found at key string: ");
                // for (const auto key : current_indices) std::cout << key << " ";
                // printf("\n");
                // auto z = 0;
                // for(const auto& key_value : current_indices){
                //     target_partial_value_counter[non_zero_indices[z]][key_value]++;
                //     z++;
                // }
                count++;
            }
        }
    }

    // printf("Round 5 Subkey: \n");
    // for(const auto& set : target_partial_value_counter){
    //     auto max_occurences = std::max(set.begin(), set.end());
    //     auto probable_key = std::distance(set.begin(), max_occurences); //get the index of the maximum probability (which is the output difference)
    //     printf("%d ", probable_key);
    // }
    // printf("\n");
    printf("Key matches: %d", count);
}

//GPT
int vector_to_hex(const std::vector<int>& vec) {
    int result = 0;
    int size = vec.size();
    
    for (int i = 0; i < size; ++i) {
        result += vec[i] * std::pow(16, size - 1 - i);
    }
    
    return result;
}

std::vector<int> hex_to_vector(int num) {
    std::vector<int> vec;

    if(num == 0){
        vec.push_back(0);
        return vec;
    }

    while (num > 0) {
        int remainder = num % 16;
        vec.insert(vec.begin(), remainder);
        num /= 16;
    }

    return vec;
}

std::vector<size_t> getRandomIndices(size_t numIndices, size_t totalIndices) {
    std::vector<size_t> result;
    std::unordered_set<size_t> indices;

    if (numIndices > totalIndices) {
        numIndices = totalIndices; // Ensure we do not request more indices than possible
    }

    std::srand(static_cast<unsigned int>(std::time(nullptr))); // Seed the random number generator

    while (indices.size() < numIndices) {
        size_t randomIndex = std::rand() % totalIndices;
        if (indices.find(randomIndex) == indices.end()) {
            indices.insert(randomIndex);
            result.push_back(randomIndex);
        }
    }

    return result;
}

//Isolate right pairs from a given set of pairs
//Techinically not able to do this but need to, to reduce compute time for demo
// std::pair<std::vector<std::pair<int, int>>, float> isolate_right_pairs(const std::vector<std::pair<int, int>>& input_pairs, const std::array<int, 16>& input_diff) {
//     std::vector<std::pair<int, int>> right_pairs;
//     auto desired_output_diff = binary_to_dec_sixteen(expected_output_diff(input_diff).first); //compute the output difference as a decimal
//     for(const auto& pair : input_pairs) {   //run the SPN for each pair and check if the output difference matches the expected output difference
//         auto ciphertext1 = spn_three_rnd(dec_to_binary_sixteen(pair.first));
//         auto ciphertext2 = spn_three_rnd(dec_to_binary_sixteen(pair.second));
//         auto output_diff = binary_to_dec_sixteen(array_xor_sixteen(ciphertext1, ciphertext2));
//         if(output_diff == desired_output_diff) right_pairs.push_back(pair);
//     }
//     // std::cout << "Number of right pairs: " << right_pairs.size() << "\n";
//     // std::cout << "Computed Output Probability: " << (double)right_pairs.size()/(double)input_pairs.size()*100 << "%\n";
//     return std::make_pair(right_pairs, float(right_pairs.size())/float(input_pairs.size())); //return the right pairs and the probability
// }

