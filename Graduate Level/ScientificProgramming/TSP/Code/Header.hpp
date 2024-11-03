//
//  Header.h
//  Heur_SFP
//
//  Created by Vincent on 11/6/15.
//
//

#ifndef Header_h
#define Header_h


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
#include <cfloat>
#include <functional>
#include <queue>

using namespace std;

#define INT_INFINITY (std::numeric_limits<int>::max())
#define DOUBLE_INFINITY (std::numeric_limits<double>::max())

#define MYRAND() (rand())
#define RAND01() (((double)(rand()%10001))/(double)10000)
#define RAND(MIN,MAX) ((MIN)+(rand()%((MAX)-(MIN))))

enum Status {Feasible, UnFeasible, OutOfRange, CycleDetected, Unknown};

struct customcomp {
private:
    double *_w;
public:
    customcomp(double *w) {_w=w;}
    bool operator() (int i,int j) { return (_w[i]+0.0001<_w[j]);}
    double operator() (int i) { return _w[i];}
    
};


template <class T, class Alloc= allocator<T> >
T RouletteWheelSelection(Alloc setofitem, double* weigth){
    assert(weigth!=NULL);
    assert(!setofitem.empty());
    auto n=setofitem.size();
    
    double sum=0.0;
    for (int i = 0; i < n; i++) sum += weigth[i];
    
    double r=sum*RAND01();
    
    T val=*setofitem.begin();
    sum=0.0;
    int k=0;
    for(auto v:setofitem){
        sum+=weigth[k];
        if(sum>=r){
            val=v;
            break;
        }
        k++;
    }
    
    
    return val;
}

template <class ForwardIterator>
ForwardIterator RandomSelect (ForwardIterator first, ForwardIterator last){
    
    if(first==last) return last;
    
    auto n = distance(first, last);
    
    advance(first, rand()%n);
    return first;
    
}


#endif /* Header_h */
