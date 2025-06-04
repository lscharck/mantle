#include <iostream>
#include <cassert>

bool smoke_test();
bool meta_test();

int main(){

    bool t1 = smoke_test();
    bool t2 = meta_test();

    assert(t1 == true);
    assert(t2 == true);

    std::cout << "Meta Test: Passed" << std::endl;

}

