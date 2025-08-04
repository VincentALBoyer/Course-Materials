//
//  Header.h
//  Heur_SFP
//
//  Created by Vincent on 11/6/15.
//
//

#pragma once


#include <ilconcert/iloenv.h>
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include <string.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <set>
#include <list>
#include <iomanip>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <limits.h>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <regex>
#include <unordered_set>
#include <unordered_map>
#include <cfloat>
#include <functional>
#include <queue>
#include <stack>
#include <iterator>

using namespace std;

#define INT_INFINITY (std::numeric_limits<int>::max())
#define DOUBLE_INFINITY (std::numeric_limits<double>::max())

#define MYRAND() (rand())
#define RAND01() (((double)(rand()%10001))/(double)10000)
#define RAND(MIN,MAX) ((MIN)+(rand()%((MAX)-(MIN))))
#define EPSILON 1e-4

enum Status {Feasible, UnFeasible, OutOfRange, CycleDetected, Unknown};

class ORUtils {
public:
    
	static bool isInteger(double x) { return fabs(x - round(x)) < EPSILON; }
	static bool isFractional(double x) { return !isInteger(x); }
	static bool isZero(double x) { return fabs(x) < EPSILON; }

	static double Floor(double x) { return std::floor(x); }
    static double FractionalPart(double x) { if (isInteger(x)) return 0.0; else return x - Floor(x); }
	static double Fractionality(double x) { return abs(0.5 - x); }

	static double coordtogeo(double val);

    // Trim from the start (left)
    static std::string ltrim(std::string s);

    // Trim from the end (right)
    static std::string rtrim(std::string s);

    // Trim from both ends
    static std::string trim(std::string s) {  return ltrim(rtrim(s));}

    static std::string getData(ifstream& input, string keyword);

    static bool gotoSection(ifstream& input, string keyword);

    static bool gotoKeyword(ifstream& input, string keyword);

    static std::string extractfilename(std::string filepath, bool withext=false);

    static std::string extractpath(std::string filepath);

    static double dotproduct(const std::vector<double>& a, const std::vector<double>& b);

};

class ORRandom {
    std::random_device rd;
    std::mt19937 gen;

    static ORRandom *_instance;

    ORRandom(unsigned int seed) : gen(seed) {
		gen = std::mt19937(rd()); // Use random_device for seeding
        // Seed the random number generator with a unique value
        gen.seed(seed);
	}


public:

    
    static ORRandom& get() {
        if (!_instance) {
            auto chrono_seed = std::chrono::system_clock::now().time_since_epoch().count();
            unsigned int seed = static_cast<unsigned int>(static_cast<uint64_t>(chrono_seed));
            _instance = new ORRandom(seed);
        }
        return *_instance;
	}


    // Delete copy/move to enforce singleton
    ORRandom(const ORRandom&) = delete;
    ORRandom& operator=(const ORRandom&) = delete;
    ORRandom(ORRandom&&) = delete;
    ORRandom& operator=(ORRandom&&) = delete;


    // Static access to the engine
    static std::mt19937& engine() {
        return get().gen;
    }

    // Uniform integer in [min, max]
    static int randint(int min, int max) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(engine());
    }

    // Uniform real in [min, max)
    static double randreal(double min = 0.0, double max = 1.0) {
        std::uniform_real_distribution<double> dist(min, max);
        return dist(engine());
    }

    // Roulette wheel selection: weights must be non-negative, not all zero
    template<typename Container>
    static size_t roulette(const Container& weights) {
        double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        if (sum <= 0.0) return 0;
        std::uniform_real_distribution<double> dist(0.0, sum);
        double r = dist(engine());
        double acc = 0.0;
        size_t idx = 0;
        for (auto w : weights) {
            acc += w;
            if (r <= acc + 0.0001) return idx;
            ++idx;
        }
        return weights.size() - 1;
    }

    // Randomly select an iterator from a container
    template<typename Container>
    static auto random_element(const Container& c) -> decltype(std::begin(c)) {
        auto size = std::distance(std::begin(c), std::end(c));
        auto it = std::begin(c);
        if (size == 0) return it;
        std::advance(it, randint(0, size - 1));
        return it;
    }

    static void destroy() {
        if (_instance) {
            delete _instance;
            _instance = nullptr;
        }
    }

};