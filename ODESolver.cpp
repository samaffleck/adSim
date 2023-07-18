#include <iostream>
#include <vector>
#include "Bed.h"

using std::cout; 
using std::endl;

int main(int argc, char **argv)
{

    Bed* myBed = new Bed();
    myBed->run_simulation();

    delete myBed;
    return 0;
}
