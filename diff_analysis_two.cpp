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

//---------------------------Core Functiions------------------------------//
void key_gen(void);
int s_box_transform(int);
int permutation_transform(int);
std::array<int, 16> run_spn(const std::array<int, 16> &); // run full SPN
void subkey_five_attack(void);
std::tuple<int,int, float> attack_partial_one(void);
std::tuple<int,int, float> attack_partial_two(void);

//-----------------Probability Functions and Generators-------------------//
std::array<std::array<int, 16>, 16> build_characteristic(void);
std::vector<std::pair<int, int>> generate_plaintext_pairs(int);
std::vector<size_t> getRandomIndices(size_t, size_t);
std::pair<int, float> expected_output_diff(int);
std::array<std::tuple<int, int, float>, 5> compute_best_differences(void);

//--------------------------Helper Functions------------------------------//
int binary_to_dec_sixteen(const std::array<int, 16> &);
int concatenate_ints(const std::array<int, 4> &inputs);
std::array<int, 16> dec_to_binary_sixteen(int);
std::array<int, 4> split_ints(int);
int vector_to_hex(const std::vector<int> &);
std::vector<int> hex_to_vector(int);

// S-Box Definition
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
    {0xF, 0x7}};

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
    {0x7, 0xF}};

// Permutation Definition
const std::map<int, int> permutation = {
    {0, 0},
    {1, 4},
    {2, 8},
    {3, 12},
    {4, 1},
    {5, 5},
    {6, 9},
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

std::array<int, 5> keys;

std::array<std::array<int, 16>, 16> prob_characteristic;

int main()
{
    key_gen();
    prob_characteristic = build_characteristic();
    auto best_diffs = compute_best_differences();
    auto [k2, k4, p1] = attack_partial_one(); //Get K2 and K4
    auto [k1, k3, p2] = attack_partial_two(); //Get K1 and K3
    printf("The fifth round subkey is: [ %d %d %d %d ] with probability %.3f%\n", k1, k2, k3, k4, p1*p2*100);
    return 0;
}

// random number generator from CPP reference
void key_gen(void)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, 1);

    for (auto &key : keys)
    {
        std::array<int, 16> key_bits;
        for (auto &bit : key_bits)
        {
            bit = distrib(gen);
        }
        key = binary_to_dec_sixteen(key_bits);
    }
    std::cout << "Keys generated: \n";
    int i = 1;
    for (const auto &key : keys)
    {
        std::cout << "Subkey " << i++ << ": ";
        std::cout << key << " = [ ";
        auto vec = split_ints(key);
        for (const auto &ele : vec)
            std::cout << ele << " ";
        std::cout << "]\n\n";
    }
}

// Take an input (that should be 4 bits) and return the output based on the S-box mapping
int s_box_transform(int input)
{
    return s_box.at(input);
}

// Take an input (that should be 16 bits) and return the output based on bit permutation mapping
int permutation_transform(int input)
{
    std::array<int, 16> input_binary = dec_to_binary_sixteen(input);
    std::array<int, 16> output;
    for (auto i = 0; i < 16; i++)
    {
        output[permutation.at(i)] = input_binary[i];
    }
    return binary_to_dec_sixteen(output);
}

// Run the SPN for a given plaintext
int run_spn(int plaintext)
{
    int output = plaintext;
    for (auto i = 0; i < 3; i++)
    {
        output ^= keys[i];               // add round key
        auto split = split_ints(output); // split the input into an array of 4 4-bit ints
        for (auto &digit : split)
            digit = s_box_transform(digit);     // apply s-box to each
        output = concatenate_ints(split);       // recombine into 16-bit int
        output = permutation_transform(output); // apply permutation transform
    }
    output ^= keys[3];               // add fourth round key
    auto split = split_ints(output); // split in 4
    for (auto &digit : split)
        digit = s_box_transform(digit); // apply s-box
    output = concatenate_ints(split);   // recombine
    output ^= keys[4];                  // add final round key
    return output;
}

// We will attack the cipher using the most probable input to output difference of:
//  0xB0B0 -> 0x8080 most probable
//  0x0B00 -> 0x0606 from the document

// Returns Partial K2 and K4 using 0x0B00 input diff
std::tuple<int,int, float> attack_partial_one(void)
{
    auto plaintext_pairs = generate_plaintext_pairs(0x0b00);              // Generate plaintext pairs
    auto random_indices = getRandomIndices(5000, plaintext_pairs.size()); // Select 5000 random indices
    std::vector<std::pair<int, int>> ciphertext_pairs;                    // Define the vector storing the output ciphertext pairs
    std::array<int, 256> subkey_counts = {0};                             // Keep count for the 256 possible partial subkeys

    for (const auto &index : random_indices)
    { // Generate 5000 ciphertext pairs for each difference
        auto ciphertext_one = run_spn(plaintext_pairs[index].first);
        auto ciphertext_two = run_spn(plaintext_pairs[index].second);
        ciphertext_pairs.push_back(std::make_pair(ciphertext_one, ciphertext_two)); // Push the ciphertext pair to the vector after running SPN on both
    }

    for (const auto &pair : ciphertext_pairs)
    {
        auto split_one = split_ints(pair.first); // Split the ciphertexts into 4 4 bit ints
        auto split_two = split_ints(pair.second);

        std::vector<int> v1 = {split_one[1], split_one[3]}; // We are only getting data from the second and fourth partial subkey, so we skip index 0 and index 2
        std::vector<int> v2 = {split_two[1], split_two[3]};

        for (auto i = 0; i < 256; i++)
        {
            auto xor1 = hex_to_vector(vector_to_hex(v1) ^ i); // Turn the two element vector into a hex value, then xor it with the partial subkey
            auto xor2 = hex_to_vector(vector_to_hex(v2) ^ i);

            for (auto &val : xor1)
                val = inverse_s_box.at(val); // Take both elements in the vector and run it through the inverse s-box
            for (auto &val : xor2)
                val = inverse_s_box.at(val);

            auto in_diff = vector_to_hex(xor1) ^ vector_to_hex(xor2); // Take the xor of both elements to calculate their last round input difference

            if (in_diff == 0x66)
            { // Check if the differential characteristic occurs
                auto current_key_array = split_ints(i);
                auto current_key = current_key_array[2] * 16 + current_key_array[3];
                subkey_counts[current_key]++; // Increment the count at whatever the current key is
            }
        }
    }

    auto max_element = std::max_element(subkey_counts.begin(), subkey_counts.end()); // Get an iterator to the largest element in the array of counts
    auto subkey = std::distance(subkey_counts.begin(), max_element);                 // The index of that iterator is the distance from the beginning
    auto partials = hex_to_vector(subkey);                                           // Convert the value to a vector so we have X -> [K2 K4]
    auto probability = float(*max_element) / 5000;                                    // Calculate the probability of the differential characteristic
    printf("Calculating from input difference of 0x0B00 (document example)\n");
    printf("Round 5 Partial K2: %d, K4: %d\n", partials[0], partials[1]);
    printf("With probability: %.3f%\n\n", probability * 100);
    return std::make_tuple(partials[0], partials[1], probability);
}

// Returns partial K1 and K3 using 0xB0B0 input difference
std::tuple<int, int, float> attack_partial_two(void)
{
    auto plaintext_pairs = generate_plaintext_pairs(0xB0B0);

    auto random_indices = getRandomIndices(5000, plaintext_pairs.size());

    std::vector<std::pair<int, int>> ciphertext_pairs;
    std::array<int, 256> subkey_counts = {0};

    for (const auto &index : random_indices)
    {
        auto ciphertext_one = run_spn(plaintext_pairs[index].first);
        auto ciphertext_two = run_spn(plaintext_pairs[index].second);
        ciphertext_pairs.push_back(std::make_pair(ciphertext_one, ciphertext_two));
    }

    for (const auto &pair : ciphertext_pairs)
    {
        auto split1 = split_ints(pair.first);
        auto split2 = split_ints(pair.second);

        std::vector<int> v1 = {split1[0], split1[2]};
        std::vector<int> v2 = {split2[0], split2[2]};

        for (auto i = 0; i < 256; i++)
        {
            auto xor1 = hex_to_vector(vector_to_hex(v1) ^ i);
            auto xor2 = hex_to_vector(vector_to_hex(v2) ^ i);

            for (auto &val : xor1)
                val = inverse_s_box.at(val);
            for (auto &val : xor2)
                val = inverse_s_box.at(val);

            auto in_diff = vector_to_hex(xor1) ^ vector_to_hex(xor2);

            if (in_diff == 0x88)
            {
                auto current_key_array = split_ints(i);
                auto current_key = current_key_array[2] * 16 + current_key_array[3];
                subkey_counts[current_key]++;
            } // Differential characteristic occurs
        }
    }

    auto max_element = std::max_element(subkey_counts.begin(), subkey_counts.end());
    auto subkey = std::distance(subkey_counts.begin(), max_element);
    auto partials = hex_to_vector(subkey);
    auto probability = float(*max_element) / 5000;
    printf("Calculataing from input difference of 0xB0B0 (from top 5)\n");
    printf("Round 5 Partial K1: %d, K3: %d\n", partials[0], partials[1]);
    printf("With probability: %.3f%\n\n", probability * 100);
    std::array<int, 2> out = {partials[0], partials[1]};
    return std::make_tuple(partials[0], partials[1], probability);
}

// Build difference characteristics for a given s-box
std::array<std::array<int, 16>, 16> build_characteristic(void)
{
    std::array<std::array<int, 16>, 16> characteristic = {0};
    for (auto i = 0; i < 16; i++)
    {
        for (auto j = 0; j < 16; j++)
        {
            auto input_diff = j ^ i;       // input difference delta_x
            auto y_1 = s_box_transform(i); // output y'
            auto y_2 = s_box_transform(j); // output y''
            auto output_diff = y_1 ^ y_2;  // get delta_y
            characteristic[input_diff][output_diff]++;
        }
    }
    printf("Differential characteristic generated:\n\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\n");
    for(auto i = 0; i<16; i++){
        printf("%d\t", i);
        for(auto j=0; j<16; j++){
            printf("%d\t", characteristic[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return characteristic;
}

// Generate a dataset of plaintext pairs for which the pairs have a given input difference
std::vector<std::pair<int, int>> generate_plaintext_pairs(int input_diff)
{
    std::vector<std::pair<int, int>> input_pairs;
    for (auto i = 0; i < 65536; i++)
    { // iterate through all possible plaintext pairs
        for (auto j = 0; j < 65536; j++)
        {
            auto current_diff = i ^ j;
            if (current_diff != input_diff)
                continue; // skip pairs that do not have the input difference
            input_pairs.push_back(std::make_pair(i, j));
        }
    }
    return input_pairs;
}

std::pair<int, float> expected_output_diff(int input_diff)
{
    auto input_array = split_ints(input_diff); // Split the input into 4 4-bit ints
    float probability = 1;                     // Initialize the probability to 1
    for (auto i = 0; i < 3; i++)
    { // Three rounds of S-boxes
        for (auto &value : input_array)
        {
            if (value == 0)
                continue;                                                                                              // Skip any zeros (unaffected S-boxes)
            auto max_prob_it = std::max_element(prob_characteristic[value].begin(), prob_characteristic[value].end()); // get the position of the maximum probability for the input difference
            probability *= float(*max_prob_it) / 16;                                                                   // multiply the probability of the output difference
            value = std::distance(prob_characteristic[value].begin(), max_prob_it);                                    // get the index of the maximum probability (which is the output difference)
        }
        auto output_diff = concatenate_ints(input_array); // Concatenate the expected output from 4 arrays into 1 int
        output_diff = permutation_transform(output_diff); // Apply permutation transform
        input_array = split_ints(output_diff);            // Split it again, and then go into next round
    }
    return std::make_pair(concatenate_ints(input_array), probability); // Return a pair containing the output difference, and the probability
}

std::array<std::tuple<int, int, float>, 5> compute_best_differences(void)
{
    std::array<std::tuple<int, int, float>, 5> best_differences;
    std::array<std::tuple<int, int, float>, 65536> all_differences;

    for (auto i = 0; i < 65536; i++)
    {                                                               // Every possible input difference
        auto [output_diff, prob] = expected_output_diff(i);         // Compute the expected output difference for the current index
        all_differences[i] = std::make_tuple(i, output_diff, prob); // Push the index, output difference, and probability to the array
    }

    std::sort(all_differences.begin(), all_differences.end(), [](std::tuple<int, int, float> &left, std::tuple<int, int, float> &right)
              {
                  return std::get<2>(left) > std::get<2>(right); // Sort the array in descending order (highest probability is the first element)
              });
    for (auto i = 1; i < 6; i++)
    { // Start at 1 because index 0 has a probability of 1
        best_differences[i - 1] = all_differences[i];
    }

    printf("Top 5 Difference Characteristics:\n");
    for (const auto& diff : best_differences){
        auto [in, out, prob] = diff;
        printf("Input Difference: %d\t\tOutput Difference: %d\t\tProbability: %.2f%\n", in, out, prob*100);
    }
    printf("\n");
    return best_differences;
}

// Convert a 16-bit binary array to a decimal number
int binary_to_dec_sixteen(const std::array<int, 16> &binary)
{
    int decimal = 0;
    for (auto i = 0; i < 16; i++)
    {
        decimal += binary[i] * (1 << (15 - i));
    }
    return decimal;
}

std::array<int, 16> dec_to_binary_sixteen(int decimal)
{
    std::array<int, 16> binary;
    for (auto i = 15; i >= 0; i--)
    {
        binary[i] = decimal % 2;
        decimal = decimal / 2;
    }
    // std::cout << "Binary value: ";
    // for (const auto& bit : binary) std::cout << bit << " ";
    return binary;
}

int concatenate_ints(const std::array<int, 4> &inputs)
{
    int result = 0;
    for (auto i = 0; i < 4; i++)
    {
        result += inputs[i] * pow(16, (3 - i));
    }
    return result;
}

std::array<int, 4> split_ints(int input)
{
    std::array<int, 4> output;
    for (int i = 0; i < 4; ++i)
        output[3 - i] = (input >> (i * 4)) & 0xF;
    return output;
}

std::vector<size_t> getRandomIndices(size_t numIndices, size_t totalIndices)
{
    std::vector<size_t> result;
    std::unordered_set<size_t> indices;

    if (numIndices > totalIndices)
    {
        numIndices = totalIndices; // Ensure we do not request more indices than possible
    }

    std::srand(static_cast<unsigned int>(std::time(nullptr))); // Seed the random number generator

    while (indices.size() < numIndices)
    {
        size_t randomIndex = std::rand() % totalIndices;
        if (indices.find(randomIndex) == indices.end())
        {
            indices.insert(randomIndex);
            result.push_back(randomIndex);
        }
    }

    return result;
}

// GPT
int vector_to_hex(const std::vector<int> &vec)
{
    int result = 0;
    int size = vec.size();

    for (int i = 0; i < size; ++i)
    {
        result += vec[i] * std::pow(16, size - 1 - i);
    }

    return result;
}

std::vector<int> hex_to_vector(int num)
{
    std::vector<int> vec;

    if (num == 0)
    {
        vec.push_back(0);
        return vec;
    }

    while (num > 0)
    {
        int remainder = num % 16;
        vec.insert(vec.begin(), remainder);
        num /= 16;
    }

    return vec;
}