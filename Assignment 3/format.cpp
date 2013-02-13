// 
// Performs various reporting functions
//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem

int main(){
     int niter = 10000;
     double t= 1000.0;
     double mx =   0.953885, l2norm=   .469796;
     printf("iteration= %d, t=% g\n",niter,t);
     printf("Max norm: %13.6e, L2norm: %13.6e\n",mx,l2norm);
     cout <<          setw(6);
     cout.setf(ios::fixed);
     cout << "iteration= " << niter << ", t= ";
     cout.unsetf(ios::fixed);
     cout.setf(ios::scientific);
     cout.precision(6);
     cout <<  t << endl;
     cout << "Max norm: " << mx << ", L2norm: " << l2norm << endl;
     return 0;
}
