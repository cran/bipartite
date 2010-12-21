// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// fitHRG - fits a hierarchical random graph (hrg) model to data
// Copyright (C) 2005-2009 Aaron Clauset
// Copyright (C) 2010 Rouven Strauss
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// See http://www.gnu.org/licenses/gpl.txt for more details.
// 
// ****************************************************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 26 October 2005
// Modified     : many, many times
//		: 27 December 2007 (cleaned up for public consumption)
//
// Modified by Rouven Strauss:
//		March  4, 2010:		modified methods for usage with bipartite networks
//		March 18, 2010:		functionality for exiting program after maxconverge time steps in which best likelihood didn't increase
//					added input parameter handling for quantitative networks
//					modified functions for usage with quantitative networks
//		April 12, 2010:		enhanced program by holding a copy of current best dendrogram instead of always writing to file immediately
//		April 13, 2010:		added further logging information
//		June 3-4, 2010:		cleaning up
//
//
// ****************************************************************************************************
// *** PROGRAM USAGE NOTES ****************************************************************************
// 
// The input to the algorithm must be a text file containing an edge list with corresponding edge weights for the graph in question; 
// vertices and edge weights are indexed by non-negative integers only, values are separated by a tab and each line is terminated by a
// carriage return. Multi-edges may appear, but will be stripped out automatically. In case of one-mode (uni-partite) networks self-loops are allowed, but have to be indicated by the flag -selfLoops.
// The graph may not be comprised of disconnected components.
//
// For instance, here is a pair of triangles linked by a single edge where all edge weights are 1 (this corresponds to a binary network):
//
// 1	2	1
// 1	3	1
// 2	3	1
// 4	5	1
// 4	6	1
// 5	6	1
// 1	4	1
//
// If the input .pairs file is formatted incorrectly, the program will crash.
//
// ****************************************************************************************************

#include <R.h>
//#include <Rmath.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"
#include "time.h"

#include "dendro.h"
#include "graph.h"
#include "rbtree.h"

using namespace std;

extern "C" { // wrapper for R

// ******** Function Prototypes ***************************************************************************

bool	markovChainMonteCarlo();
string	num2str(const unsigned int);
bool	parseCommandLine(int argc, char * argv[]);
bool	readInputFile();
void	recordModules();
void	recordNamesLUT();

// ******** Structures and Constants **********************************************************************

struct ioparameters {
	int		n_a;			// number of A vertices in input graph
	int		n_b;			// number of B vertices in input graph
	int		n;			// total number of vertices in input graph (i.e. n_a + n_b)
	int		m;			// number of edges in input graph
	unsigned int	maxconverge;		// maximum number of steps without increase of best likelihood
	double		temperature;		// SA temperature for computing modules
	string		d_dir;			// working directory
	string		f_in;			// name of input file (*.pairs)
	string		f_simA;			// name of input file containing similarity information about A vertices
	string		f_simB;			// name of input file containing similarity information about B vertices
	string		f_dg;			// name of output file
	string		f_dg_info;		// name of output information file
	string		f_ordA;			// name of output order file for A vertices
	string		f_ordB;			// name of output order file for B vertices
	string		f_modules;		// name of output file for modules
	string		f_nestedModules;	// name of output file for nested modules
	string		f_stat;			// name of output statistics file
	string		f_pairs;		// name of output random graph file
	string		f_namesLUT;		// name of output names LUT file
	string		s_scratch;		// filename sans extension
	string		s_tag;			// user defined filename tag
	string		start_time;		// time simulation was started
	int		timer;			// timer for reading input
	bool		flag_filename;		// flag indicating whether -filename invoked
	bool		flag_steps;		// flag indicating whether -steps invoked
	bool		flag_timer;		// flag for when timer fires
	bool		flag_compact;		// compact the Lxy file
	bool		flag_bipartite;		// flag indicating whether -bipartite invoked
	bool		flag_selfLoops;		// flag indicating whether -selfLoops invoked
};

// ******** Global Variables ******************************************************************************

ioparameters	ioparm;				// program parameters
rbtree		namesLUT;			// look-up table; translates input file vertex names to graph indices
rbtree		reverseNamesLUT;		// reverse look-up table; translates graph indices to input file vertex names
dendro*		d;				// hrg data structure
dendro*		bestDendro;			// dendrogram with current best likelihood
unsigned int	t;				// number of time steps max = 2^32 ~ 4,000,000,000
unsigned int	t_bestDendro;			// number of time steps until best dendrogram was found
double		temperature;			// current SA temperature
double		dTemperature;			// SA temperature step
double		averageStartTemperature;	// average temperature at which increase of
int		averageDenominator;		// auxiliary variable
unsigned int	converge;			// current number of steps without increase of best likelihood
double		minTemperature;			// minimum SA temperature
short int	billionCount;			// counts number of billion steps
double		bestL;				// best likelihood found so far
int		nrOfRecordBreakings;		// counts number of maximum found
unsigned int	period  = 10000;		// number of MCMC moves to do before writing stuff out; default: 10000
double*		Likeli;				// holds last k hrg likelihoods
char*		task;				// task to execute (find modules)	
MTRand		mtr;				// Mersenne Twister random number generator instance

// ******** Main Loop *************************************************************************************

int identifyModules(int* r_argc, char* argv[]) {

	ioparm.n_a			= 0;				//
	ioparm.n_b			= 0;				//
	ioparm.n			= 0;				//
	ioparm.temperature		= 1e-5;				//
	ioparm.timer   			= 20;				//
	ioparm.flag_filename		= false;			//
	ioparm.flag_steps		= false;			//
	ioparm.flag_compact	 	= true;				//
	ioparm.flag_bipartite		= false;			//
	ioparm.flag_selfLoops		= false;			//
	ioparm.s_tag   			= "";				//
	ioparm.maxconverge		= 0;				//
	minTemperature			= 1e-999;			//
	string input   			= "";				//
	t				= 1;				//
	billionCount			= 0;				//
	nrOfRecordBreakings		= 0;				//
	task				= "computeModules";		//
	time_t t1			= time(&t1);			//
	int argc			= *r_argc;			//

	if (parseCommandLine(argc, argv)) {
		d = new dendro(task);							// make the dendro-graph structure for computing
		ioparm.start_time = asctime(localtime(&t1));
		Likeli = new double [period];						// allocate space

		if (!readInputFile()) {
			 cout << "!! Error: Malformed input file.\n"; return 0;
		}

 		bestDendro		= d->deepCopy();				// make the dendro-graph structure holding a copy of dendrogram with best likelihood
 		bestL			= d->getLikelihood();				// store current likelihood
 
		temperature		= ioparm.temperature;
		dTemperature		= (temperature - minTemperature) / (double)(ioparm.maxconverge);

		cout << "\nstep\tL\tbest L\tMC step\t\tdL\ttemperature\n";
		while(converge < ioparm.maxconverge) {						// likelihood did increase during the last maxconverge steps
			if (!(markovChainMonteCarlo())) { return 0; }
			if (t >= 1000000000) { billionCount++; t = 0; }			// rollover step count
		}

		recordModules();							// write dendrogram to file

		delete bestDendro;
		bestDendro = NULL;
		delete d;
		d = NULL;

		return 1;
	}
	else {
		return 0;
	} 

}

// ******** Function Definitions **************************************************************************

bool markovChainMonteCarlo() {

	double  dL;
	bool    flag_taken;

	// Because moves in the dendrogram space are chosen (Monte Carlo) so that we sample dendrograms 
	// with probability proportional to their likelihood, a likelihood-proportional sampling of 
	// the dendrogram models would be equivalent to a uniform sampling of the walk itself. We would
	// still have to decide how often to sample the walk (at most once every n steps is recommended)
	// but for simplicity, the code here simply runs the MCMC itself. To actually compute something
	// over the set of sampled dendrogram models (in a Bayesian model averaging sense), you'll need
	// to code that yourself.
	
	// do 'period' MCMC moves before doing anything else
	for (unsigned int i=0; i<period; i++) {
		
		if (!(d->monteCarloMove(dL, flag_taken, temperature, bestL))) {				// make a MCMC move
			cout << "!! Error: failed to make monte carlo move" << endl;
			return false;
		}

		Likeli[i] = d->getLikelihood();								// get likelihood of this D given G
		
		if (Likeli[i] > bestL) {
			
			bestL = Likeli[i];								// store the current best likelihood
			
			if(averageDenominator == 0) {
				averageStartTemperature	= temperature;
				averageDenominator	= 1;
			}
			else {
				averageStartTemperature *= averageDenominator;
				averageStartTemperature += temperature;
				averageDenominator++;
				averageStartTemperature /= (double)(averageDenominator);
			}

			temperature	= averageStartTemperature + mtr.randExc()*(ioparm.temperature - averageStartTemperature);
			dTemperature	= (temperature - minTemperature) / (double)(ioparm.maxconverge);
			
			nrOfRecordBreakings++;								// increment count of record-breakings

			delete bestDendro;
			bestDendro = d->deepCopy();							// copy d to bestDendro
			t_bestDendro = t;
			
			converge = 0;								// reset convergence counter
		}
		else { 
			if(temperature - dTemperature >= minTemperature) temperature -= dTemperature;

			converge++;									// increase convergence counter
		}

	
		// Write some stuff to standard-out to describe the current state of things.
		if ((ioparm.flag_compact and ((i+(t-period)) % (8192) == 1)) or !ioparm.flag_compact) {
			cout << "[" << billionCount << " | " << t << "]\t" << Likeli[i];
			cout << "   \t(" << bestL << ")\t";
			if (flag_taken) { cout << "*\t"; } else { cout << " \t"; }
			cout << "\t" << dL << " \t" << temperature << endl;
		}

		t++;
	}

	d->refreshLikelihood();										// corrects floating-point errors O(n)

	return true;
}

// ********************************************************************************************************

string num2str(const unsigned int input) {
	// input must be a positive integer
	unsigned int temp = input;
	string str  = "";
	if (input == 0) {
		str = "0";
	}
	else {
		while (temp != 0) {
			str  = char(int(temp % 10)+48) + str;
			temp = (unsigned int)temp/10;
		}
	}
	return str;
}

// ********************************************************************************************************

bool parseCommandLine(int argc, char * argv[]) {
	int argct = 1;
	string temp, ext;
	string::size_type pos;

	if (argc==1) {
		cout << "\n  -- Hierarchical Module Identification --\n";
		cout << "  by Rouven Strauss (copyright 2010-2020)\n\n";
		cout << "  based on the algorithm \n";
		cout << "\n  -- Hierarchical Random Graphs --\n";
		cout << "  by Aaron Clauset (copyright 2005-2008)\n\n";
		cout << "  Flags:\n";
		cout << "  -filename <file>               (required) input .pairs graph file\n";
		cout << "  -steps <integer>               (required) maximum number of steps without\n";
		cout << "                                            increase of best likelihood (for module detection)\n";
		cout << "  -bipartite                     (required if input graph is bipartite)\n";
		cout << "  -selfLoops                     (required if self-loops shall be allowed)\n";
		cout << "  -label <string>                    (optional) label for this run\n";
		cout << "\n";
		cout << "  examples:\n";
		cout << "  ./identifyModules -filename graph.pairs -steps 1000000\n";
		cout << "  ./identifyModules -filename graph.pairs -steps 1000000 -label test -selfLoops\n";
		cout << "  ./identifyModules -filename graph.pairs -steps 1000000 -bipartite\n";
		cout << "\n";
		return false;
		
	} else {
		while (argct < argc) {
			temp = argv[argct];
			
			if (temp == "-label") {
				argct++;
				ioparm.s_tag = argv[argct];
			}
			else if (temp == "-compact") {
				ioparm.flag_compact = true;
			}
			else if (temp == "-filename") {
				ioparm.flag_filename = true;
				argct++;
				temp = argv[argct];
				ext = ".pairs";
				pos = temp.find(ext,0);
				if (pos == string::npos) {
					cout << "!! Error: Input file must claim to be .pairs format.\n";
					return false;
				}
				ioparm.f_in = temp;
				ext = "/";
				pos = string::npos;
				for (unsigned int i=0; i < temp.size(); i++) {
					if (temp[i] == '/') {
						pos = i;
					}
				}
				if (pos != string::npos) {
					ioparm.d_dir = temp.substr(0, pos+1);
					temp = temp.substr(pos+1,temp.size()-pos-1);
				}
				// now grab the filename sans extension for building outputs files
				for (unsigned int i=0; i < temp.size(); i++) {
					if (temp[i] == '.') {
						pos = i;
					}
				}
				ioparm.s_scratch = temp.substr(0,pos);
				
			}
			else if (temp == "-bipartite") {
				ioparm.flag_bipartite = true;
			}
			else if (temp == "-selfLoops") {
				ioparm.flag_selfLoops	= true;
			}
			else if (temp == "-steps") {
				ioparm.flag_steps = true;
				argct++;
				if(atoi(argv[argct]) < 0) {
					cout << "!! Error: -steps argument has to be >= 0!" << endl;
					return false;
				}
				else {									// hand over an integer representing the maximum number
					ioparm.maxconverge = atoi(argv[argct]);				// of steps without increase of best likelihood before exiting
				}
			}
			else if (temp == "-experimental") {
				task = "computeModulesExperimental";
			}
			else {
				cout << "Warning: ignored argument " << argct << " : " << temp << endl;
			}
			argct++;
		}
	}

	ioparm.f_namesLUT = ioparm.d_dir + ioparm.s_scratch + "-names.lut";
	
	if (ioparm.s_tag != "")    {
		ioparm.s_scratch += "_" + ioparm.s_tag;
	}

	if (ioparm.flag_filename)         {
		ioparm.f_stat     = ioparm.d_dir + ioparm.s_scratch + "-L.xy";
	}

	if(ioparm.flag_bipartite && ioparm.flag_selfLoops) {
		cout << "!! Error: there are no bipartite graphs with self-loops!" << endl;
		cout << "          Please adjust flags!" << endl;
		return false;
	}
	else if(ioparm.flag_bipartite) {
		cout << ">> Graph is bipartite" << endl;
	}
	else if(ioparm.flag_selfLoops) {
		cout << ">> Graph is unipartite and self loops are allowed." << endl;
	}
	else {
		cout << ">> Graph is unipartite and self loops are forbidden." << endl;
	}

	if(!ioparm.flag_filename) {
		cout << "!! Error: flag -filename required!" << endl;
		return false;
	}

	if(!ioparm.flag_steps) {
		cout << "!! Error: -steps has to be invoked with appropriate parameters!" << endl;
		return false;
	}

	return true;
}

// ********************************************************************************************************

bool readInputFile() {

	int n_a, n_b, n, m, vertex_i, vertex_j, edgeWeight, virtualVertex_i, virtualVertex_j, countBVertices, sumEdgeWeight;
	n_a = n_b = n = m = countBVertices = sumEdgeWeight = 0;
	elementrb *item;
	time_t t1; t1 = time(&t1);
	time_t t2; t2 = time(&t2);

	// First, we scan through the input file to create a list of unique vertex names
	// (which we store in the namesLUT), and a count of the number of edges.
	cout << "\n>> input file scan ( " << ioparm.f_in << " )" << endl;

	cout << ">> starting to count edges" << endl;
	cout << ">> edges: [0]"<<endl;

	ifstream fscan0(ioparm.f_in.c_str(), ios::in);						// read input

	while (fscan0 >> vertex_i >> vertex_j >> edgeWeight) {					// read friendship pair (vertex_i,vertex_j)
		if (vertex_i == vertex_j && !ioparm.flag_selfLoops) {
			if(ioparm.flag_bipartite) {
				cout << "!! Error: there are no bipartite graphs with self-loops!" << endl;
			}
			else {
				cout << "!! Error: self-loops are forbidden!" << endl;
				cout << "!!        Activate -selfLoops if self-loops are desired." << endl;
			}
			return false;
		}
		else {
			m++;									// count number of edges
			sumEdgeWeight += edgeWeight;						// compute total sum of edge weights
			if (namesLUT.findItem(vertex_i) == NULL) {
				namesLUT.insertItem(vertex_i, n_a);
				n_a++;								// increment amount of A vertices
			}
			if(!ioparm.flag_bipartite) {
				if (namesLUT.findItem(vertex_j) == NULL) {
					namesLUT.insertItem(vertex_j, n_a);
					n_a++;								// increment amount of A vertices
				}
			}
		}
			
		if (t2-t1>ioparm.timer) {							// check timer; if necessary, display
			cout << ">> edges: ["<<m<<"]"<<endl;
			t1 = t2; ioparm.flag_timer = true;					//
		}										//
		t2=time(&t2);									//
	}

	fscan0.close();
	
	if(ioparm.flag_bipartite) {

		countBVertices = n_a;
	
		ifstream fscan1(ioparm.f_in.c_str(), ios::in);

		while (fscan1 >> vertex_i >> vertex_j >> edgeWeight) {					// read friendship pair (vertex_i,vertex_j)
			if (vertex_i != vertex_j) {
				if (namesLUT.findItem(vertex_j) == NULL) {
					namesLUT.insertItem(vertex_j, countBVertices);
					countBVertices++;
					n_b++;								// increment amount of B vertices
				}
			}
			else {
				cout << "!! Error: there are no bipartite graphs with self-loops." << endl;
				return false;
			}
		}

		fscan1.close();

	}
	
	cout << ">> total amount of edges: ["<<m<<"]"<<endl;
	cout << ">> sum of edge weights: ["<<sumEdgeWeight<<"]"<<endl;

	d->g = new graph(n_a, n_b, ioparm.flag_bipartite);						// make new graph with (n_a + n_b) vertices

	// Finally, we reparse the file and add edges to the graph
	m			= 0;
	ioparm.flag_timer	= false;								// reset timer
	bool aToB		= true;

	ifstream fin(ioparm.f_in.c_str(), ios::in);

	while (fin >> vertex_i >> vertex_j >> edgeWeight) {
		m++;
		if (!ioparm.flag_bipartite || vertex_i != vertex_j) {
			item = namesLUT.findItem(vertex_i); virtualVertex_i = item->value;
			item = namesLUT.findItem(vertex_j); virtualVertex_j = item->value;
			if (!(d->g->doesLinkExist(virtualVertex_i, virtualVertex_j))) {
				if (!(d->g->addLink(virtualVertex_i, virtualVertex_j, edgeWeight, aToB))) {
					cout << "!! Error: couldn't insert edge (" << vertex_i << " " << vertex_j << " " << edgeWeight <<")" << endl;
					return false;
				}
				else if (d->g->getName(virtualVertex_i) == "") {
					d->g->setName(virtualVertex_i, num2str(vertex_i));
				}
			}
			if (!(d->g->doesLinkExist(virtualVertex_j, virtualVertex_i))) {
				if (!(d->g->addLink(virtualVertex_j, virtualVertex_i, edgeWeight, !aToB))) {
					cout << "!! Error: couldn't insert edge (" << vertex_i << " " << vertex_j << " " << edgeWeight <<")" << endl;
					return false;
				}
				else if (d->g->getName(virtualVertex_j) == "") {
					d->g->setName(virtualVertex_j, num2str(vertex_j));
				}
			}
		}
	}

	fin.close();

	// computes the expected edge weights for all vertex pairs
	d->g->computeExpectedEdgeWeight();

	ioparm.m	= d->g->getNumLinks();							// store number of directional edges created
	ioparm.n_a	= d->g->getNumAVertices();							// store number of A vertices used
	ioparm.n_b	= d->g->getNumBVertices();							// store number of B vertices used
	ioparm.n	= d->g->getNumVertices();							// store total number of vertices used
	cout << ">> A vertices: ["<<ioparm.n_a<<"]"<<endl;
	cout << ">> B vertices: ["<<ioparm.n_b<<"]"<<endl;
	cout << ">> total amount of vertices: ["<<ioparm.n<<"]" << endl;
	
	recordNamesLUT();									// record names LUT to file for future reference
	
	d->buildDendrogram(ioparm.flag_bipartite);

	return true;
}

// ********************************************************************************************************

void recordModules() {

	time_t t1;

	// write to files
	ioparm.f_dg = ioparm.d_dir + ioparm.s_scratch + ".hrg";
	ioparm.f_ordA = ioparm.d_dir + ioparm.s_scratch + ".ordA";
	ioparm.f_ordB = ioparm.d_dir + ioparm.s_scratch + ".ordB";
	ioparm.f_modules = ioparm.d_dir + ioparm.s_scratch + ".mod";
	ioparm.f_nestedModules = ioparm.d_dir + ioparm.s_scratch + ".mod_exp";

	//bestDendro->recordDendrogramStructure(ioparm.f_dg);
	if(!d->recordOrderAndModules(reverseNamesLUT, ioparm.f_ordA, ioparm.f_ordB, ioparm.f_modules, ioparm.f_nestedModules)) {
		cout << "!! Error: failed to record order and module files" << endl;
		return;
	}
	if(!strcmp(task, "computeModulesExperimental")) {
		remove(ioparm.f_modules.c_str());
		rename(ioparm.f_nestedModules.c_str(), ioparm.f_modules.c_str());
	}
	else if(!strcmp(task, "computeModules")) {
		remove(ioparm.f_nestedModules.c_str());
	}
	
	// write statistics about hrg to file
	ioparm.f_dg_info = ioparm.d_dir + ioparm.s_scratch + ".info";

	t1 = time(&t1); 
	ofstream fout(ioparm.f_dg_info.c_str(), ios::trunc);
	fout << "---HIERARCHICAL-RANDOM-GRAPH---\n";
	fout << "StartTime	: " << ioparm.start_time;
	fout << "InputFile	: " << ioparm.f_in    << "\n";
	fout << "Directory	: " << ioparm.d_dir   << "\n";
	fout << "---Basic-Graph-Information---\n";
	if(ioparm.flag_bipartite) {
		fout << "Bipartite	: yes\n";
	}
	else {
		fout << "Bipartite	: no\n";
	}
	fout << "A vertices	: " << ioparm.n_a << "\n";
	fout << "B vertices	: " << ioparm.n_b << "\n";
	fout << "Total amount of vertices	: " << ioparm.n << "\n";
	fout << "Edges		: " << ioparm.m/2 << "\n";
	fout << "---HRG-Information---\n";
	fout << "OutputTime	: " << asctime(localtime(&t1));
	fout << "NumStepsMCMC	: " << t << "\n";
	fout << "NumPrevBests	: " << nrOfRecordBreakings-1 << "\n";
	fout << "Likelihood	: " << bestL << "\n";
	fout << "HRG		: " << ioparm.s_scratch + ".hrg" << "\n";
	fout << "InfoFile	: " << ioparm.s_scratch + ".info"  << "\n";
	fout.close();
	
	cout << ">> recorded dendrogram info to file: " << ioparm.f_dg_info << endl;
	
	return;
}

// ********************************************************************************************************

void recordNamesLUT() {
	keyValuePair *head, *prev;
	
	head = namesLUT.returnTreeAsList();
	while (head != NULL) {
		reverseNamesLUT.insertItem(head->y, head->x);
		prev = head;
		head = head->next;
		delete prev;
	}
	head = NULL; prev = NULL;
	
	elementrb *item;
	ofstream fout(ioparm.f_namesLUT.c_str(), ios::trunc);
	fout << "virtual\treal\n";
	for (int i=0; i<ioparm.n; i++) {
		item = reverseNamesLUT.findItem(i);
		fout << i << "\t" << item->value << "\n";
	}
	fout.close();
	
	cout << ">> recorded names look-up table to file: " << ioparm.f_namesLUT << endl;
	
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************

} // end wrapper for R
