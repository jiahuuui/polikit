// main.cpp

#include <iostream>
#include "nf.h"

using namespace std;
using namespace polikit;

void box::cvalue(){
        lx=10.0;
        ly=82.8;
        lz=3.46;
        natom=90;

    }

int main(){

    class box t;
    t.helpmessage();
    t.cvalue();
    // class box t;
    // t.helpmessage();

    cout<<t.lx<<t.ly<<t.natom<<endl;
    return 0;
}
