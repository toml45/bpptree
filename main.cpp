//
// Created by tomtom on 5/14/2024.
//

#include <iostream>
#include <new>
#include "ADS_set.h"

int main() {
    //std::vector<int> v{1,2,3};
    ADS_set<int> *test = new ADS_set<int>{11,15,14,12,13,17,18,20,22,25,26,27,28,29,33,1};
    test->dump();
    //std::cout << test->count(2) << std::endl;
    delete test;
}