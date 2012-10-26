// Process command line arguments
// 
//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include "types.h"
using namespace std;

void cmdLine(int argc, char *argv[], int& n, int& BSZ, int& do_blas, int& do_hand_blocked, int& do_unblocked, int &nreps){
/// Command line arguments
 // Default value of the domain sizes
 static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"b", required_argument, 0, 'b'},
        {"nreps", required_argument, 0, 'r'},
        {"do_blas", no_argument, 0, 's'},
        {"do_hand_blocked", no_argument, 0, 'h'},
        {"do_unblocked", no_argument, 0, 'u'},
 };

 // Set default values
    do_blas=0, do_hand_blocked=0, do_unblocked=0;
    n=16, nreps=5, BSZ = 16;

    // Process command line arguments
 int ac;
 for(ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"n:b:r:shu",long_options,NULL)) != -1){
        switch (c) {

	    // Size of the matrix
            case 'n':
                n = atoi(optarg);
                break;

	    // Blocking factor
            case 'b':
                BSZ = atoi(optarg);
                break;

	    // Use BLAS in blocked variant
            case 's':
                do_blas = 1;
                break;

	    // Use hand coded blocking variant
            case 'h':
                do_hand_blocked = 1;
                break;

	    // Use hand coded blocking variant
            case 'u':
                do_unblocked = 1;
                break;

	    // # of repititions
            case 'r':
                nreps = atoi(optarg);
		if (nreps < 1){
		    cerr << "nreps must be > 0. Exiting....\n\n";
		    exit(-1);
		}
                break;

	    // Error
            default:
                cout << "Usage: mmpy [-n <matrix dim>] [-b <cache blocking factor]  [-r <# repititions>]";
                cout << "\n\t    ";
                cout << " [-s <use BLAS>] [-h <use hand blocking>] [-u <unblocked>]\n\t";
                cout << endl;
                exit(-1);
            }
    }
 }
}
