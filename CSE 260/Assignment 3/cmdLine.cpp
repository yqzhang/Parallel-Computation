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

void cmdLine(int argc, char *argv[], int& n, int& stats_freq, int& plot_freq, int& px, int& py, int &noComm, int &niter){
/// Command line arguments
 // Default value of the domain sizes
 static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"stats-freq", required_argument, 0, 's'},
        {"plot", required_argument, 0, 'p'},
	{"px", required_argument, 0, 'x'},
	{"py", required_argument, 0, 'y'},
	{"niter", required_argument, 0, 'i'},
	{"nocomm", no_argument, 0, 'k'},

 };
    // Process command line arguments
 int ac;
 for(ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"n:x:y:i:j:s:kp:",long_options,NULL)) != -1){
        switch (c) {

	    // Size of the computational box
            case 'n':
                n = atoi(optarg);
                break;
	   //
           // X processor geometry
            case 'x':
                px = atoi(optarg);
                break;

            // X processor geometry
            case 'y':
                py = atoi(optarg);
                break;

            // # of iterations
	    // Use this option control the number of mesh sweeps
            case 'i':
                niter = atoi(optarg);
                break;


	    // Print various statistics
            case 's':
                stats_freq = atoi(optarg);
                break;


	    // Plot the excitation variable
            case 'p':
                plot_freq = atoi(optarg);
                break;
	    //
            // Shut off communication
            case 'k':
                noComm = 1;
                break;

	    // Error
            default:
                cout << "Usage: apf [-n <domain size>] [-i <# iterations>]";
                cout << "\n\t    ";
                cout << " [-s <stats frequency>[-p <plot frequency>]\n\t";
		cout << "     [-tx <x processor geometry> [-ty <y processor geometry] [-k <no communication>]" << endl;
                cout << endl;
                exit(-1);
            }
    }
 }
}
