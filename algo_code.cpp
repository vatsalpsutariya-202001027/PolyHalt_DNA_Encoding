#include<bits/stdc++.h>
#include <iostream>
#include <fstream>
#include<math.h>
#include <random>
// #include <algorithm>
using namespace std;
#define endl "<br>"
std:: ifstream fin("input.txt");
std:: ofstream fout("output.txt");
// Declare variables for input parameters
int N, X, Y, H;
int k = 0;
int precautionCount = 0;
int precautionLimit = 1e5;
int totalCount = 0;
vector<string> all;
vector<string> filtered;
string curr;
string nuc = "ACGT";
unordered_map<char, char> reverseComplement;
unordered_map<string, int> prohibited;

std::random_device rd;
std::mt19937 gen(rd());

int a = 1;
int b = 1001;
std::uniform_int_distribution<int> distribution(a, b);
int randomNum1 = distribution(gen);
int randomNum2 = distribution(gen);


// Initialize the reverse complement map for DNA bases
void initReverseComplement() {
  reverseComplement['A'] = 'T';
  reverseComplement['T'] = 'A';
  reverseComplement['G'] = 'C';
  reverseComplement['C'] = 'G';
  return;
}

// Function to get the reverse complement of a DNA sequence
string getReverseComplement(string &s) {
  string res;
  for (auto a : s) {
    res.push_back(reverseComplement[a]);
  }
  return res;
}
std::string generateRandomDNA(int N) {
    // srand(time(0));

    // DNA alphabet
    const std::string DNA = "ATCG";

    // Initialize the DNA sequence with random characters
    std::string randomDNA;
    for (int i = 0; i < N; ++i) {
        int randomIndex = distribution(gen) % DNA.length();
        randomDNA += DNA[randomIndex];
    }

    return randomDNA;
}

// Cost function to measure the "uniqueness" of a DNA sequence
// Example: the more unique characters, the lower the cost
double costFunction(const std::string& dna) {
    int uniqueCharacters = 0;
    for (char base : "ATCG") {
        if (dna.find(base) != std::string::npos) {
            uniqueCharacters++;
        }
    }
    return static_cast<double>(uniqueCharacters);
}

// Stochastic Local Search Algorithm to generate unique DNA sequence
std::string stochasticLocalSearch(int N, int maxIterations) {
    // srand(time(0));

    std::string currentDNA = generateRandomDNA(N);

    for (int iteration = 1; iteration <= maxIterations; ++iteration) {
        // Evaluate the cost of the current DNA sequence
        double currentCost = costFunction(currentDNA);

        // Generate a random neighbor DNA sequence
        std::string neighborDNA = generateRandomDNA(N);

        // Evaluate the cost of the neighbor DNA sequence
        double neighborCost = costFunction(neighborDNA);

        // Accept the neighbor DNA sequence if it has a lower cost
        if (neighborCost < currentCost) {
            currentDNA = neighborDNA;
        }
    }

    return currentDNA;
}

void generateRandom(){
  srand(time(0));
  for(int i = 0 ; i < 100000 ; i ++){
    all.push_back(stochasticLocalSearch(N, distribution(gen)%15));
  }
  return;
}

// Function to generate all possible combinations of DNA sequences of length N
void generateCombinations() {
  // if(totalCount >= 1e5) return;
  if (curr.size() == N) {
    // totalCount++;
    all.emplace_back(curr); // Store the generated sequence in the 'all' vector
    return;
  }
  if(totalCount >= 1e5) return;
  // srand((int)time(0));
  int st = distribution(gen);
  // fout  << st << endl;
  st %= 4;
  for (int i = 0 ; i < 4 ; i ++) {
    char n = nuc[(i+st)%4];
    curr.push_back(n); // Add a DNA base to the current sequence
    generateCombinations(); // Recursively generate the next base
    curr.pop_back(); // Backtrack and remove the last added base
  }
  return;
}

// Function to check for repeating patterns within a DNA sequence
bool check(int degree, string &s) {
  int n = s.size();
  for (int i = 0; i < n; i++) {
    int startOne = i;
    int startTwo = degree + i;
    string temp1 = s.substr(startOne, degree); // Extract a substring of length 'degree'
    string temp2 = "";
    if (startTwo < n) {
      temp2 = s.substr(startTwo, degree); // Extract another substring starting 'degree' bases ahead
    }
    if (temp1 == temp2) return false; // Check if the two substrings are equal, indicating a repeating pattern
  }
  return true;
}

// Function to check if a DNA sequence is homopolymer-free up to a certain degree
bool homopolymerFree(string &s, int limit) {
  for (int degree = 1; degree <= limit; degree++) {
    if (!check(degree, s)) return false; // Check for repeating patterns up to the given 'limit'
  }
  return true;
}

// Function to check if a DNA sequence is free of secondary structures
bool secondaryStructureFree(string &s, int threshold) {
  string rc = s; // Create a reverse complement of the input DNA sequence
  reverse(rc.begin(), rc.end());
  for (auto &c : rc) {
    c = reverseComplement[c]; // Replace each base with its complement in the reverse sequence
  }
  int store = N;
  N = s.size();
  int maxi = -1;
  vector<vector<int>> dp(N + 1, vector<int>(N + 1, 0)); // Initialize a 2D array for dynamic programming
  for (int i = N; i >= 0; i--) {
    for (int j = N; j >= 0; j--) {
      if (i == N || j == N) {
        dp[i][j] = 0;
        continue;
      }
      if (s[i] == rc[j]) {
        dp[i][j] = 1 + dp[i + 1][j + 1]; // If the bases match, increment the score
        maxi = max(maxi, dp[i][j]); // Track the maximum score
      } else {
        dp[i][j] = 0; // If the bases don't match, reset the score to 0
      }
    }
  }
  N = store;
  if (maxi > threshold) {
    return false; // If the maximum score is above the threshold, the sequence has secondary structures
  }
  return true;
}

// Function to check if a DNA sequence has balanced GC content
bool balancedGC(string &s) {
  unordered_map<char, int> mp;
  for (auto a : s) {
    mp[a]++; // Count the occurrences of each base
  }
  if(mp['A'] + mp['T'] == mp['G'] + mp['C']) return true; // Check if A-T and G-C base pairs are balanced
  return false;
}

// Function to calculate the Edit distance between two DNA sequences
int findEditDis(string &s, string &t) {
  
  // Code here
  int n = s.size();
  int m = t.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
  
    for (int i = 0; i <= n; i++) {
      for (int j = 0; j <= m; j++) {
        if (i == 0 || j == 0) {
          dp[i][j] = i + j;
        } else if (s[i - 1] == t[j - 1]) {
          dp[i][j] = dp[i - 1][j - 1];
        } else {
          dp[i][j] = min(
            min(
            1 + dp[i - 1][j],  // Insert.
            1 + dp[i][j - 1]),  // Remove.
            1 + dp[i - 1][j - 1]);  // Replace.
        }
      }
    }
  
    return dp[n][m];
}

bool valid = false;

// Function to find a set of DNA codewords with a minimum Edit distance
void getXCodewordsPrint(vector<string> &allcurr, int i, vector<string> &curr) {
  if (valid == true) return; // If a valid set is found, stop searching
  if(precautionCount >= precautionLimit){
    return ;
  }
  if (i == 0) {
    for (auto a : curr) {
      for (auto b : curr) {
        precautionCount++;
        if (a != b && findEditDis(a, b) < H) return; // Check Edit distance for each pair in the current set
      }
    }
    valid = true; // A valid set with minimum Edit distance is found
    for (auto a : curr) fout  << a << endl; // Print the set of codewords
    fout  << endl << endl;
    return;
  }
  if (valid == true) return;
  for (auto a : allcurr) {
    int flag = 0;
    for (auto b : curr) {
      precautionCount++;
      if (findEditDis(a, b) < H) {
        flag = 1;
        break;
      }
    }
    if (flag) {
      continue; // Skip sequences that don't meet the Edit distance requirement
    }
    curr.push_back(a); // Add a sequence to the current set
    getXCodewordsPrint(allcurr, i - 1, curr); // Recursively search for more codewords
    curr.pop_back(); // Backtrack by removing the last added sequence
  }
  return;
}

void getXCodewords(vector<string> &allcurr, int i, vector<string> &curr) {
  if (valid == true) return; // If a valid set is found, stop searching
  if(precautionCount >= precautionLimit){
    return ;
  }
  if (i == 0) {
    for (auto a : curr) {
      for (auto b : curr) {
        precautionCount++;
        if (a != b && findEditDis(a, b) < H) return; // Check Edit distance for each pair in the current set
      }
    }
    valid = true; // A valid set with minimum Edit distance is found
    // for (auto a : curr) fout  << a << endl; // Print the set of codewords
    // fout  << "<br><br>";
    return;
  }
  if (valid == true) return;
  for (auto a : allcurr) {
    int flag = 0;
    for (auto b : curr) {
      precautionCount++;
      if (findEditDis(a, b) < H) {
        flag = 1;
        break;
      }
    }
    if (flag) {
      continue; // Skip sequences that don't meet the Edit distance requirement
    }
    curr.push_back(a); // Add a sequence to the current set
    getXCodewords(allcurr, i - 1, curr); // Recursively search for more codewords
    curr.pop_back(); // Backtrack by removing the last added sequence
  }
  return;
}

int main() {
  srand(time(0));
    // Mersenne Twister 19937 generator

    // Define a uniform distribution for integers between a and b (inclusive)
    
    // Generate two random numbers

  // fout  << "Enter the size of codeword (N): ";
  fin >> N;
  // fout  << "read";
  // N = 6; X = 3; Y = 2; H = 3;
  // fout  << "Enter X such that you want to prevent homopolymers up to X degrees: ";
  fin >> X;
  
  // fout  << "Enter the maximum length allowed for the patterns for secondary structures (Y): ";
  fin >> Y;
  // fout  << "Enter the minimum Edit / Edit distance required (H): ";
  fin >> H;
  initReverseComplement();
  totalCount = 0;
  if(N <= 10) generateCombinations();
  else{
    generateRandom();
  }
  // fout  << all.size();
  // return 0;
  for (auto &a : all) {
    if (homopolymerFree(a, X) && balancedGC(a)) {
      if (secondaryStructureFree(a, Y)) {
        prohibited[getReverseComplement(a)]++; // Store the reverse complement of the sequence
        if (prohibited[a] == 0) {
          filtered.push_back(a);
          // fout << a << endl;
        }
      }
    }
  }
  // shuffling the final codeword vector to obtain most random results
  int siz = filtered.size();
  for (int i = 0; i < siz - 1; i++) {
    int j = i + rand() % (siz - i);
    swap(filtered[i], filtered[j]);
  }
  // for(auto a: all) fout  << a << endl;
  // return 0;
  vector<string> curr;
  int low = 1;
  int high = 300; // Find a set of codewords with a minimum Edit distance
  
    bool exists = false;
  while(low <= high){
    int mid = (low + high) / 2;
    // fout  << "Checking " << mid;
    valid = false;
    precautionCount = 0;
    getXCodewords(filtered, mid+1, curr);
    curr.clear();
    if(valid){
      exists = true;
      low = mid + 1;
      // fout  << " - Positive<br>"; 
    }
    else {
      high = mid-1;
      // fout  << " - Negative<br>";
    }
  }
  
  int mid = (low + high) / 2;
  fout << "<br>";
  if(exists) fout  << "Total codewords are : " << (low + high) / 2 << endl;
  fout << "<br><br><br>";
  valid = false;

  precautionCount = 0;
  getXCodewordsPrint(filtered, mid, curr);
  if (!valid) {
    fout << -1;
    // fout  << "Sorry, no codewords can be formed with these requirements.";
  }
  else{
    
    fout << endl << "You are having information density of: " << endl;
    double ans = log((low + high)/2)/log(2);
    ans /= N;
    fout  << ans << endl;
  }
  return 0;
}
