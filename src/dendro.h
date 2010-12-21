// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// dendro.h - hierarchical random graph (hrg) data structure
// Copyright (C) 2005-2009 Aaron Clauset
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
// Created      : fitHRG:	26 October 2005 - 7 December 2005
//		: consensusHRG:	19 April 2006
//		: predictHRG:	21 June 2006
//
// Modified     : fitHRG:	23 December 2007 (cleaned up for public consumption)
//		: consensusHRG:	19 May 2007
//		: 		19 May 2008 (cleaned up for public consumption)
//		: predictHRG:	23 December 2007 (cleaned up for public consumption)
//
// Modified by Rouven Strauss:
//		March 8-9, 2010:	modified algorithm for usage with bipartite networks
//		March 10-11, 2010:	modified likelihood function for usage with bipartite networks
//		March 13, 2010:		merged dendro.h files of fitHRG, consensusHRG and predictHRG into one
//		March 15, 2010:		changed computing of likelihood:
//						likelihood is now computed by summing up the internal vertices' contribution to likelihood which is defined as:
//							-(n-1)			if Lna*Rnb+Lnb*Rna = 0
//							e/(Lna*Rnb+Lnb*Rna)	else
//		March 16, 2010:		added methods for order file containing:
//						- row and column ordering information
//		March 18, 2010:		modified methods for usage with weighted networks (invocation of addLink() in makeRandomGraph(), quant)
//		March 22-23, 2010:	added margin totals approach
//		March 24, 2010:		integrated different likelihood approaches used until now
//					cleaned up code
//					modified method for order file output in order to include also information about modules (position, size and level of depth)
//		March 25-26, 2010:	(visweb) added functionality for plotting modules
//		March 31, 2010:		preparatory work for modularity computation
//		April 6, 2010:		changed likelihood to ew/nL_nR and debugged expected margin likelihood
//		April 8, 2010:		added methods for modularity computation
//		April 13, 2010:		added method for deep copying dendrogram
//					addded method for merging modules
//		April 14, 2010:		finished methods for modularity computation
//		May 10, 2010:		changed approach for computing best modules (outside/inside module approach)
//		June 3-4, 2010:		cleaning up
//
// ****************************************************************************************************
// 
// Maximum likelihood dendrogram data structure. This is the heart of the HRG algorithm: all
// manipulations are done here and all data is stored here. The data structure uses the separate
// graph data structure to store the basic adjacency information (in a dangerously mutable way).
// 
// ****************************************************************************************************

#if !defined(dendro_INCLUDED)
#define dendro_INCLUDED

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <float.h>

#include "MersenneTwister.h"
#include "graph.h"
#include "rbtree.h"

using namespace std;

// ********************************************************************************************************
// ******** Basic Structures ******************************************************************************

#if !defined(list_INCLUDED)
#define list_INCLUDED
class list {
public:
	int	x;						// stored elementd in linked-list
	list*	next;						// pointer to next elementd
	list(); ~list();
};
list::list()  { x = -1; next = NULL; }
list::~list() {}
#endif

enum {DENDRO, GRAPH, LEFT, RIGHT, PARTITION_A, PARTITION_B};
#if !defined(block_INCLUDED)
#define block_INCLUDED
struct block		{ double x; int y; };
#endif
struct ipair		{ int x; int y; short int t; string sp; };
struct edgeCountTriple	{ int e; int e_w; double e_w_expect; };

// ********************************************************************************************************
// ******** Internal Edge Class ***************************************************************************

// The usefulness of this data structure is to provide an easy to way maintain the set of internal edges,
// and the corresponding splits, in the dendrogram D. It allows for the selection of a random internal 
// edge in O(1) time, and it takes O(1) time to update its structure given an internal move. 
// This structure does not provide any means to directly manipulate the splits, but does allow them to be replaced.
// A split has the form "int.int...int#int.int...int", where all ints on the left side of the # are in the left
// partition and all ints on the right side of the # marker are in the right partition defined by the split.

class interns {
private: 
	ipair*	edgelist;		// list of internal edges represented
	int**	indexLUT;		// table of indices of internal edges in edgelist
	int	nrOfLegalRandomEdges;	//
	int*	legalRandomEdges;	//
	int	q;			// number of internal edges
	int	count;			// (for adding edges) edgelist index of new edge to add
	MTRand	mtr;			// Mersenne Twister random number generator instance

public:
	interns(const int);
	~interns();
	
	bool	addEdge(const int, const int, const short int);						// add an internal edge, O(1)
	ipair*	getEdge(const int);									// returns the ith edge of edgelist, O(1)
	ipair*	getRandomEdge();									// returns a uniformly random internal edge, O(1)
	bool	swapEdges(const int, const int, const short int, const int, const int, const short int);// swaps two edges, O(1)
	void	printEdgeList();									// writes edgelist to terminal
};

// ********************************************************************************************************

interns::interns(const int n)  {
	q         = n;
	count     = 0;
	edgelist  = new ipair  [q];
	indexLUT  = new int*   [q+1];
	for (int i=0; i<(q+1); i++) {
		indexLUT[i]    = new int [2];
		indexLUT[i][0] = indexLUT[i][1] = -1;
	}
}
interns::~interns() {
	delete [] edgelist;
	for (int i=0; i<(q+1); i++) { delete [] indexLUT[i]; }
	delete [] indexLUT;
}

// ********************************************************************************************************

bool interns::addEdge(const int new_x, const int new_y, const short int new_type) {
	// This method adds a new edge (i,j,t,sp) to the list of internal edges. After checking that the inputs
	// fall in the appropriate range of values, it records the new edgelist index in the indexLUT and then
	// puts the input values into that edgelist location.
	if(count < q and new_x >= 0 and new_x < (q+1) and new_y >= 0 and new_y < (q+2) and (new_type == LEFT or new_type == RIGHT)) {
		if(new_type == LEFT) {
			indexLUT[new_x][0] = count;
		}
		else {
			indexLUT[new_x][1] = count;
		}
		edgelist[count].x = new_x;
		edgelist[count].y = new_y;
		edgelist[count].t = new_type;
		count++;
		return true;
	}
	else {
		return false;
	}
}

// ********************************************************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getEdge(const int i) {
	return &edgelist[i];
}

// ********************************************************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getRandomEdge() {

	int nr = (int)(floor((double)(q)*mtr.randExc()));

	return &edgelist[nr];
}

// ********************************************************************************************************

bool interns::swapEdges(const int one_x, const int one_y, const short int one_type, const int two_x, const int two_y, const short int two_type) {
	// The moves on the dendrogram always swap edges, either of which (or both, or neither) can by internal
	// edges. So, this method mirrors that operation for the internal edgelist and indexLUT.
	
	int index, jndex, temp;
	bool one_isInternal = false;
	bool two_isInternal = false;
	
	if(one_x >= 0 and one_x < (q+1) and two_x >= 0 and two_x < (q+1) and (two_type == LEFT or two_type == RIGHT) and 
	    one_y >= 0 and one_y < (q+2) and two_y >= 0 and two_y < (q+2) and (one_type == LEFT or one_type == RIGHT)) {
		
		if(one_type              == LEFT)	{ temp = 0; } else { temp = 1; }
		if(indexLUT[one_x][temp] >  -1  )	{ one_isInternal = true;       }
		if(two_type              == LEFT)	{ temp = 0; } else { temp = 1; }
		if(indexLUT[two_x][temp] >  -1  )	{ two_isInternal = true;       }
		
		if(one_isInternal and two_isInternal) {
			if(one_type == LEFT)  { index = indexLUT[one_x][0]; } else { index = indexLUT[one_x][1]; }
			if(two_type == LEFT)  { jndex = indexLUT[two_x][0]; } else { jndex = indexLUT[two_x][1]; }
			temp			= edgelist[index].y;
			edgelist[index].y	= edgelist[jndex].y;
			edgelist[jndex].y	= temp;
			
		} else if(one_isInternal) {
			if(one_type == LEFT)	{ index = indexLUT[one_x][0]; indexLUT[one_x][0] = -1; }
			else			{ index = indexLUT[one_x][1]; indexLUT[one_x][1] = -1; }
			edgelist[index].x	= two_x;
			edgelist[index].t	= two_type;
			if(two_type == LEFT)	{ indexLUT[two_x][0] = index; } else { indexLUT[two_x][1] = index; } // add new
			
		} else if(two_isInternal) {
			if(two_type == LEFT)	{ index = indexLUT[two_x][0]; indexLUT[two_x][0] = -1; }
			else			{ index = indexLUT[two_x][1]; indexLUT[two_x][1] = -1; }
			edgelist[index].x	= one_x;
			edgelist[index].t	= one_type;
			if(one_type == LEFT)	{ indexLUT[one_x][0] = index; } else { indexLUT[one_x][1] = index; } // add new
		} else { } // else neither is internal

		return true;
	} else { return false; }
}

// ********************************************************************************************************

void interns::printEdgeList() {
	for (int i=0; i<q; i++) {
		cout << "(" << edgelist[i].x << " " << edgelist[i].y << " ";
		if(edgelist[i].t == LEFT)  {
			cout << "L) ";
		}
		else if(edgelist[i].t == RIGHT) {
			cout << "R) ";
		}
		else {
			cout << "?) ";
		}
	}
	cout << endl;
	return;
}

// ********************************************************************************************************
// ******** Tree elementd Class ***************************************************************************

class elementd {
public:
	short int	type;				// either DENDRO or GRAPH
	short int	partition;			// either PARTITION_A or PARTITION_B
	double		Lcont;				// likelihood contribution of this internal vertex
	double		p;				// probability p_i that an edge exists between L and R subtrees
	int		e;				// number of edges between L and R subtrees
	int		e_w;				// sum of edge weights between L and R subtrees
	int		e_w_total;			// sum of edge weights between and within L and R subtrees
	double		e_w_expect;			// sum of "expected" edge weights between L and R subtrees computed by margin totals
	double		e_w_expect_total;		// sum of "expected" edge weights between and within L and R subtrees computed by margin totals
	int		n_a;				// number of leafs of partition A in subtree rooted here
	int		n_b;				// number of leafs of partition B in subtree rooted here
	int		n;				// total number of leafs in subtree rooted here
	int		nrOfFurtherModules;		// indicates whether this internal vertex contains more than one module (the one formed by itself)
	int		label;				// subtree label: smallest leaf index
	int		index;				// index in containing array
	
	elementd   	*M;				// pointer to parent vertex
	elementd   	*L;				// pointer for L subtree
	elementd   	*R;				// pointer for R subtree
	
	elementd(); ~elementd();
};
elementd::elementd()  {
	type			= DENDRO;
	partition		= PARTITION_A;
	label		= index = -1;
	e		= e_w	= 0;
	Lcont		= p	= 0.0;
	e_w_expect	= 0.0;
	n	= n_a	= n_b	= 0;
	M	= L	= R	= NULL;
	nrOfFurtherModules	= -1;
}
elementd::~elementd() {}

// ********************************************************************************************************
// ******** Dendrogram Class ******************************************************************************

class dendro {
	
private:
	elementd*	root;				// root of the dendrogram
	elementd*	internal;			// array of n-1 internal vertices (the dendrogram D)
	elementd*	leaf;				// array of n   leaf vertices (the graph G)
	bool		isBipartite;			// bool indicating whether graph G is bipartite
	int		n_a;				// number of leaf vertices in partition A
	int		n_b;				// number of leaf vertices in partition B
	int		n;				// number of leaf vertices to allocate
	int		nrOfModules;			// number of modules without further nested modules
	int		internalCount;			// needed for recordDendrogramStructure
	char*		task;				// task to execute
	int		sumEdgeWeight;			// total sum of edge weights
	interns*	d;				// list of internal edges of dendrogram D
	list**		paths;				// array of path-lists from root to leaf
	double		L;				// likelihood of graph G given dendrogram D
	MTRand		mtr;				// Mersenne Twister random number generator instance
	rbtree		subtreeL, subtreeR;		// trees for computeEdgeCount() method
	
	void		binarySearchInsert(elementd*, elementd*);								// insert vertex i according to binary search property
	list*		binarySearchFind(const double);										// return path to root from leaf
	edgeCountTriple*	computeEdgeCount(const int, const short int, const int, const short int);			// compute number of edges between two internal subtrees
	elementd*	findCommonAncestor(list**, const int, const int);							// find internal vertex of D that is common ancestor of i,j
	void		printSubTree(elementd*);										// display the subtree rooted at z
	list*		reversePathToRoot(const int);										// return reverse of path to leaf from root
	void		QsortMain(block*, int, int);										// quicksort methods
	int		QsortPartition(block*, int, int, int);
	bool		setValues(int, int, int, double, interns*, dendro*);									// set values (n_a, n_b, sumEdgeWeight, root) of dendrogramm
	void		setBackNrOfFurtherModules(elementd*);								// set back number of modules of all vertices in subtree
	int		setNrOfFurtherModules(elementd*);									// set number of modules for each internal vertex
	int		getInternalVertexEdgeWeightSum(elementd*);								// computes sum(x->e_w) of an internal vertex, where x iterates over all of its children
	int		setTotalEdgeWeight(elementd*);										// sets e_w_total for each internal vertex
	double		setTotalExpectedEdgeWeight(elementd*);									// sets e_w_expect_total for each internal vertex
	double		controlNrOfModules(elementd*);										// method in order to control whether order of nrOfFurtherModules in dendrogram is valid
	void		setNrOfFurtherModules(elementd*, int, bool, bool);							// sets nrOfFurtherModules for each internal vertex (-1 if vertex within module, 0 if vertex forms module, 1 if vertex outside of module)
	double		computeLcont(elementd*);									// computes likelihood contribution of vertices (expect of y and x) involved in MCMC move
	elementd*	getCopyOfLeaves();											// return copy of 'leaf' array
	elementd*	getCopyOfInternals(elementd*);										// return copy of 'internal' array
	list*		mergeLists(list*, list*);										// merge two already ordered lists into one ordered list
	void		deleteList(list*);											// delete list
	list*		recordOrderAndModules(rbtree&, ofstream&, ofstream&, ofstream&, ofstream&, elementd*, const int, const int);	// write order files and module files
	list*		getInternalVertexIndicesWithinModules();								// returns a list with the indices of the internal vertices which are head of a module or within one

public:
	graph*		g;													// underlying G (dangerously accessible)
	
	dendro(char*); ~dendro();													// constructor / destructor
	void		buildDendrogram(const bool);										// build dendrogram from g
	double		getLikelihood();											// return likelihood of G given D
	bool		importDendrogramStructure(const string);								// read dendrogram structure from file
	dendro*		deepCopy();												// return copy of current dendrogram
	bool		monteCarloMove(double&, bool&, const double, const double);								// make single MCMC move
	void		refreshLikelihood();											// force refresh of likelihood value
	void		printDendrogram();											// write dendrogram structure to terminal
	void		makeRandomGraph();											// make random G from D
	void		recordGraphStructure(const string);									// record G structure to file
	bool		recordOrderAndModules(rbtree&, const string, const string, const string, const string);				// invokes private recordModuleOrder method
	void		recordDendrogramStructure(const string);								// invoke private recursive method recordDendrogramStructure in order to record dendrogram structure to file
	void		setTask(char*);												// set task
};

// ******** Dendrogram Methods ****************************************************************************

dendro::dendro(char* d_task) {
	root		= NULL;
	internal	= NULL;
	leaf		= NULL;
	d		= NULL;
	paths		= NULL;
	g		= NULL;
	task		= d_task;
}
dendro::~dendro() {
	list *curr, *prev;
	
	//if(g		!= NULL) { delete g;           g		= NULL; }    	// O(m)
	
	if(internal	!= NULL) { delete [] internal; internal  = NULL; }    		// O(n)
	if(leaf	!= NULL) { delete [] leaf;     leaf	     = NULL; }    	// O(n)
	if(d		!= NULL) { delete d;           d         = NULL; }    		// O(n)
	if(paths	!= NULL) {
		for (int i=0; i<n; i++) {
			curr = paths[i]; while (curr != NULL) {
				prev = curr;
				curr = curr->next;
				delete prev;
				prev = NULL;
			}
			paths[i] = NULL;
		}
		delete [] paths;
	}
	paths = NULL;
}

// ********************************************************************************************************
// *** private methods ************************************************************************************
// ********************************************************************************************************

void dendro::binarySearchInsert(elementd* x, elementd* y) {
	if(y->p < x->p) {					// go to left subtree
		if(x->L == NULL) { 				// check if left subtree is empty
			x->L = y;				// make x left child
			y->M = x;				// make y parent of child
			return;
		}
		else { binarySearchInsert(x->L, y); }
	} else {						// go to right subtree
		if(x->R == NULL) { 				// check if right subtree is empty
			x->R = y;				// make x right child
			y->M = x;				// make y parent of child
			return;
		}
		else { binarySearchInsert(x->R, y); }
	}
	return;
}

// ********************************************************************************************************

list* dendro::binarySearchFind(const double v) {
	list *head, *tail, *newlist;
	elementd *current = root;
	bool flag_stopSearch = false;
	
	while (!flag_stopSearch) {										// continue until we're finished
		newlist    = new list;										// add this vertex to the path
		newlist->x = current->label;
		if(current == root) { head       = newlist; tail = head;    }
		else                 { tail->next = newlist; tail = newlist; }
		if(v < current->p) {										// now try left subtree
			if(current->L->type == GRAPH) { flag_stopSearch = true; }
			else { current = current->L; }
		} else {											// else try right subtree
			if(current->R->type == GRAPH) { flag_stopSearch = true; }
			else { current = current->R; }
		}
	}
	return head;
}

// ********************************************************************************************************

edgeCountTriple* dendro::computeEdgeCount(const int x, const short int xtype, const int y, const short int ytype) {
	// This method computes the number of edges that cross between the subtree internal[x]
	// and the subtree internal[y] and the sum of their weights.
	// To do this, we use an array X[1..n] integers which take values -1 if X[i] is in the
	// subtree defined by internal[x], +1 if X[i] is in the subtree internal[y], and 0
	// otherwise. Taking the smaller of the two sets, we then scan over the edges attached
	// to that set of vertices and count the number of endpoints we see in the other set.

	int nX, nY;
	bool flag_go			= true;
	edgeCountTriple* count	= new edgeCountTriple;
	count->e			= 0;
	count->e_w			= 0;
	count->e_w_expect		= 0.0;
	const short int k		= 1 + DENDRO + GRAPH;

	elementd* curr;
	
	// --- First, we push the leaf vertices in the L and R subtrees into balanced binary tree
	//     structures so that we can search them quickly later on.
	if(xtype == GRAPH) {								// default case, subtree X is size 1
		subtreeL.insertItem(x,-1);						// insert single vertex as member of left subtree
		nX		 = 1;							// 
	}
	else {
		curr		= &internal[x];						// explore subtree X, O(|X|)
		curr->type	= k+1;							//
		nX		= 0;							//
		while (flag_go) {
			if(curr->index == internal[x].M->index) {
				internal[x].type= DENDRO;
				flag_go		= false;
			} else {
				if(curr->type == k+1 and 
				    curr->L->type == GRAPH) {				// - is it time, and is left child a graph vertex?
					subtreeL.insertItem(curr->L->index, -1);
					curr->type = k+2;	  			//
					nX++;						//
				}
				if(curr->type == k+2 and curr->R->type == GRAPH) {	// - is it time, and is right child a graph vertex?
					subtreeL.insertItem(curr->R->index, -1);
					curr->type = k+3;	  			//
					nX++;						//
				}
				if(curr->type == k+1) {				// - go left
					curr->type	= k+2;				//
					curr		= curr->L;			//
					curr->type	= k+1;
				}
				else if(curr->type == k+2) {				// - else go right
					curr->type = k+3;				// 
					curr       = curr->R;				// 
					curr->type = k+1;	
				}
				else {							// - else go up a level
					curr->type = DENDRO;				// 
					curr       = curr->M;				// 
					if(curr == NULL) {
						flag_go = false;
						cout << "X exit: reached null parent" << endl;
					}
				}
			}
			if(nX > n) { cout << "error! nX > n\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n" << endl; break; }
		}
	}
	
	if(ytype == GRAPH) {								// default case, subtree Y is size 1
		subtreeR.insertItem(y,1);						// insert vertex as single member of right subtree
		nY = 1;									// 
	} else {
		flag_go = true;
		curr		= &internal[y];						// explore subtree Y, O(|Y|)
		curr->type	= k+1;							//
		nY		= 0;							//
		while (flag_go) {	
			if(curr->index == internal[y].M->index) {
				internal[y].type = DENDRO;
				flag_go = false;
			}
			else {
				if(curr->type == k+1 and curr->L->type == GRAPH) {	// - is it time, and is left child a graph vertex?
					subtreeR.insertItem(curr->L->index, 1);
					curr->type = k+2;				//
					nY++;						// 
				}
				if(curr->type == k+2 and 
				    curr->R->type == GRAPH) {				// - is it time, and is right child a graph vertex?
					subtreeR.insertItem(curr->R->index, 1);
					curr->type = k+3;				// 
					nY++;						//
				}
				if(curr->type == k+1) {				// - look left
					curr->type	= k+2;				// 
					curr		= curr->L;				// 
					curr->type	= k+1; }
				else if(curr->type == k+2) {				// - look right
					curr->type	= k+3;				// 
					curr		= curr->R;			// 
					curr->type	= k+1;
				} else {						// - else go up a level
					curr->type	= DENDRO;			// 
					curr		= curr->M;			// 
					if(curr == NULL) {
						flag_go	= false;
					}
				}
			}
			if(nY > n) { cout << "error! nY > n \n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n" << endl; break; }
		}
	}

	// --- Now, we take the smaller subtree and ask how many of its emerging edges have their
	//     partner in the other subtree. If we are computing the likelihood using the margin totals
	//     time is O(|X| * |Y|), else it is O(|X| log |X|)
	edge* current;
	int*  treeListL;
	int*  treeListR;

	if(nX < nY) {									// subtreeL is smaller
		treeListL = subtreeL.returnArrayOfKeys();
		treeListR = subtreeR.returnArrayOfKeys();
		for (int i=0; i<nX; i++) {
			current = g->getNeighborList(treeListL[i]);
			while (current != NULL) {					// loop over each of its neighbors v_j
				if(subtreeR.findItem(current->x) != NULL) {
					count->e++;
					count->e_w += current->weight;
				}
				current = current->next;				// to see if v_j is in X
			}								// 
			subtreeL.deleteItem(treeListL[i]);

			for(int j = 0; j < nY; j++) {
				count->e_w_expect += g->getExpectedEdgeWeight(treeListL[i], treeListR[j]);
			}
		}
		delete [] treeListL;
		for (int i=0; i<nY; i++) { subtreeR.deleteItem(treeListR[i]); }
		delete [] treeListR;
	} else {									// subtreeR is smaller
		treeListR = subtreeR.returnArrayOfKeys();
		treeListL = subtreeL.returnArrayOfKeys();
		for (int i=0; i<nY; i++) {
			current = g->getNeighborList(treeListR[i]);
			while (current != NULL) {					// loop over each of its neighbors v_j
				if(subtreeL.findItem(current->x) != NULL) {
					count->e++;
					count->e_w += current->weight;
				}
				current = current->next;				// to see if v_j is in Y
			}								// 
			subtreeR.deleteItem(treeListR[i]);

			for(int j = 0; j < nX; j++) {
				count->e_w_expect += g->getExpectedEdgeWeight(treeListR[i], treeListL[j]);
			}
		}
		delete [] treeListR;
		for (int i=0; i<nX; i++) { subtreeL.deleteItem(treeListL[i]); }
		delete [] treeListL;
	}

	return count;
}

// ********************************************************************************************************

elementd* dendro::findCommonAncestor(list** paths, const int i, const int j) {
	list* headOne = paths[i];
	list* headTwo = paths[j];
	elementd* lastStep;
	while (headOne->x == headTwo->x) {
		lastStep = &internal[headOne->x];
		headOne  = headOne->next;
		headTwo  = headTwo->next;
		if(headOne == NULL or headTwo == NULL) { break; }
	}
	return lastStep;			// Returns address of an internal vertex; do not deallocate
}

// ********************************************************************************************************

void dendro::printSubTree(elementd *z) {
	if(z != NULL) {
		if(z->type == GRAPH) {
			cout << "[" << z->label << "]" << endl; //"\t(" << z->L << " " << z->R << ") - " << z->M <<  endl;
			return;
		} else if(z->type == DENDRO) {
			cout <<  "(p = " << z->p << "\te = " << z->e << "\tnL = " << z->L->n << "\tnR = " << z->R->n << "\tlabel = " << z->label << ")\tinternal[" << z->index << "]" << endl; // "\t(" << z->L << " " << z->R << ") - " << z->M <<  endl;
			cout << "L "; printSubTree(z->L); cout << endl;
			cout << "R "; printSubTree(z->R); cout << endl;
		} else {
			cout <<  "(p = " << z->p << "\te = " << z->e << "\tnL = " << z->L->n << "\tnR = " << z->R->n << "\tlabel = " << z->label << ")\tinternal[" << z->index << "] " << z->type << endl;
			cout << "L "; printSubTree(z->L); cout << endl;
			cout << "R "; printSubTree(z->R); cout << endl;
		}
	}
	return;
}

// ********************************************************************************************************

list* dendro::reversePathToRoot(const int leafIndex) {
	list *head, *subhead, *newlist;
	head = subhead = newlist = NULL;
	elementd *current = &leaf[leafIndex];
	
	while (current != NULL) {				// continue until we're finished
		newlist       = new list;			// add this vertex to the path
		newlist->x    = current->index;
		newlist->next = NULL;
		if(head == NULL) { head    = newlist; }
		else              { subhead = head;    head = newlist; head->next = subhead; }
		current = current->M;
	}
	return head;
}

// ********************************************************************************************************

void dendro::QsortMain (block* array, int left, int right) {
	if(right > left) {
		int pivot = left;
		int part  = QsortPartition(array, left, right, pivot);
		QsortMain(array, left,   part-1);
		QsortMain(array, part+1, right  );
	}
	return;
}

// ********************************************************************************************************

int dendro::QsortPartition (block* array, int left, int right, int index) {
	block p_value, temp;
	p_value.x = array[index].x;
	p_value.y = array[index].y;
	
	// swap(array[p_value], array[right])
	temp.x		= array[right].x;
	temp.y		= array[right].y;
	array[right].x	= array[index].x;
	array[right].y	= array[index].y;
	array[index].x	= temp.x;
	array[index].y	= temp.y;
	
	int stored	= left;
	for (int i=left; i<right; i++) {
		if(array[i].x <= p_value.x) {
			// swap(array[stored], array[i])
			temp.x		= array[i].x;
			temp.y		= array[i].y;
			array[i].x	= array[stored].x;
			array[i].y	= array[stored].y;
			array[stored].x	= temp.x;
			array[stored].y	= temp.y;
			stored++;
		}
	}
	// swap(array[right], array[stored])
	temp.x		= array[stored].x;
	temp.y		= array[stored].y;
	array[stored].x	= array[right].x;
	array[stored].y	= array[right].y;
	array[right].x	= temp.x;
	array[right].y	= temp.y;
	
	return stored;
}

// ********************************************************************************************************

bool dendro::setValues(int d_n_a, int d_n_b, int d_sumEdgeWeight, double d_L, interns* d_d, dendro* toCopy) {
	n_a		= d_n_a;			// number of vertices in Partition A
	n_b		= d_n_b;			// number of vertices in Partition B
	n		= d_n_a + d_n_b;		// total number of vertices in graph
	L		= d_L;
	sumEdgeWeight	= d_sumEdgeWeight;		// total sum of edge weights
	g		= toCopy->g;
	d		= new interns(n-2);      	// allocate memory for internal edges of D, O(n)
	
	ipair* tempPair;

	for(int i = 0; i < n-1; i++) {
		tempPair = d_d->getEdge(i);		// returns address; do not deallocate
		d->addEdge(tempPair->x, tempPair->y, tempPair->t);
	}

	leaf		= toCopy->getCopyOfLeaves();
	internal	= toCopy->getCopyOfInternals(leaf);

	root		= &internal[0];

	return true;
}

// ********************************************************************************************************
// *** following private methods needed only for fitHRG ***************************************************
// ********************************************************************************************************

void dendro::setBackNrOfFurtherModules(elementd* vertex) {

	if(vertex->type == DENDRO) {
		setBackNrOfFurtherModules(vertex->L);
		setBackNrOfFurtherModules(vertex->R);
		vertex->nrOfFurtherModules = 0;
	}
}

// ********************************************************************************************************

// Method computes for each internal vertex whether there are more modules within its subtrees
// than the one its two subtrees form together
// O(n)
int dendro::setNrOfFurtherModules(elementd* vertex) {

	if (vertex->L->type == GRAPH || vertex->R->type == GRAPH) {
		nrOfModules++;
		return 1;
	}
	else {
		vertex->nrOfFurtherModules = setNrOfFurtherModules(vertex->L) + setNrOfFurtherModules(vertex->R);
		return vertex->nrOfFurtherModules;
	}
}

// ********************************************************************************************************

// Method computes for an internal vertex the sum of all x->e_w, where x is a subtree of the internal vertex,
// e.g. getInternalVertexEdgeWeightSum(root) = sumEdgeWeight
int dendro::getInternalVertexEdgeWeightSum(elementd* vertex) {

	if(vertex->type == DENDRO) { return vertex->e_w + getInternalVertexEdgeWeightSum(vertex->L) + getInternalVertexEdgeWeightSum(vertex->R); }
	return 0;
}

// ********************************************************************************************************

int dendro::setTotalEdgeWeight(elementd* vertex) {
	if(vertex->type == DENDRO) {
		vertex->e_w_total = setTotalEdgeWeight(vertex->L) + setTotalEdgeWeight(vertex->R) + vertex->e_w;
		return vertex->e_w_total;
	}
	else return 0;
}

// ********************************************************************************************************

double dendro::setTotalExpectedEdgeWeight(elementd* vertex) {
	if(vertex->type == DENDRO) {
		vertex->e_w_expect_total = setTotalEdgeWeight(vertex->L) + setTotalEdgeWeight(vertex->R) + vertex->e_w_expect;
		return vertex->e_w_expect_total;
	}
	else return 0;
}

// ********************************************************************************************************

double dendro::controlNrOfModules(elementd* vertex) {
	if(vertex->type == GRAPH) return (double)(vertex->nrOfFurtherModules);
	else {
		double v = max(controlNrOfModules(vertex->L), controlNrOfModules(vertex->R));
		if(v > vertex->nrOfFurtherModules) { cout << "!! Error: corrupt dendrogram" << endl; }
		if(v == 0 && vertex->nrOfFurtherModules == 0) { cout << "!! Error: corrupt dendrogram" << endl; }
		return max(v, (double)vertex->nrOfFurtherModules);
	}
}

// ********************************************************************************************************

void dendro::setNrOfFurtherModules(elementd* vertex, int value, bool setLcont, bool completely) {

	if(vertex->type == DENDRO) {

		if(value == 1 && (vertex->L->type == GRAPH || vertex->R->type == GRAPH)) {
			if(setLcont) {
				if(vertex->nrOfFurtherModules == 1) {
					vertex->Lcont = -vertex->Lcont;
				}
				if(vertex->e_w == 0) {
					vertex->Lcont = (double)(-(n-1));
				}
			}
			vertex->nrOfFurtherModules = 0;
			value = -1;
		}
		else {
			if(setLcont) {
				if(value == -1 && vertex->nrOfFurtherModules == 1) {
					if(vertex->e_w == 0) {
						vertex->Lcont = (double)(-(n-1));
					}
					else {
						vertex->Lcont = -vertex->Lcont;
					}
				}
				else if(value == 1 && vertex->nrOfFurtherModules == -1) {
					if(vertex->e_w == 0) {
						vertex->Lcont = ((double)(vertex->e_w) - vertex->e_w_expect) / (double)(sumEdgeWeight);
					}
					else {
						vertex->Lcont = -vertex->Lcont;
					}
				}
			}
			vertex->nrOfFurtherModules = value;
		}

		if(!(value == -1 && vertex->L->nrOfFurtherModules == -1) || completely) {
			setNrOfFurtherModules(vertex->L, value, setLcont, completely);
		}
		if(!(value == -1 && vertex->R->nrOfFurtherModules == -1) || completely) {
			setNrOfFurtherModules(vertex->R, value, setLcont, completely);
		}
	}
}

// ********************************************************************************************************

double dendro::computeLcont(elementd* vertex) {

	if(vertex->L->type == GRAPH || vertex->R->type == GRAPH) {
		return 0.0;
	}

	double lcont_new;

	if(vertex->e_w == 0) {
		if(vertex->nrOfFurtherModules != 1) {																	// if current internal vertex with e_w = 0 is within a module
			if(vertex->Lcont != (double)(-(n-1))) cout << "!! Error: internal vertex [within module with unconnected subtrees] with wrong likelihood contribution" << endl;	// its likelihood contribution has to be -(n-1)
			lcont_new = ((double)(vertex->e_w) - vertex->e_w_expect) / (double)(sumEdgeWeight);
		}
		else {																					// if current internal vertex with e_w = 0 is outside a module
			if(vertex->L->n_a * vertex->R->n_b + vertex->L->n_b * vertex->R->n_a != 0 && vertex->Lcont == (double)(-(n-1))) {
				cout << "!! Error: internal vertex [outside module with unconnected subtrees] with wrong likelihood contribution" << endl;	// its likelihood contribution has to greater than -(n-1)
				cout << vertex->L->n_a * vertex->R->n_b + vertex->L->n_b * vertex->R->n_a << endl;
			}
			lcont_new = (double)(n-1);
		}

		return computeLcont(vertex->L) + computeLcont(vertex->R) + lcont_new + vertex->Lcont;
	}
	else {
		return computeLcont(vertex->L) + computeLcont(vertex->R) + 2 * vertex->Lcont;
	}
}

// ********************************************************************************************************

elementd* dendro::getCopyOfLeaves() {

	elementd* result = new elementd[n];

	for(int i = 0; i < n; i++) {
		result[i].type		= leaf[i].type;
		result[i].partition	= leaf[i].partition;
		result[i].n_a		= leaf[i].n_a;
		result[i].n_b		= leaf[i].n_b;
		result[i].n		= leaf[i].n;
		result[i].nrOfFurtherModules = leaf[i].nrOfFurtherModules;
		result[i].label		= leaf[i].label;
		result[i].index		= leaf[i].index;
	}

	return result;
}

// ********************************************************************************************************

elementd* dendro::getCopyOfInternals(elementd* leaves) {

	elementd* result = new elementd[n-1];

	for(int i = 0; i < n-1; i++) {
		result[i].type		= internal[i].type;
		result[i].partition	= internal[i].partition;
		result[i].Lcont		= internal[i].Lcont;
		result[i].p		= internal[i].p;
		result[i].e		= internal[i].e;
		result[i].e_w		= internal[i].e_w;
		result[i].e_w_total	= internal[i].e_w_total;
		result[i].e_w_expect	= internal[i].e_w_expect;
		result[i].e_w_expect_total = internal[i].e_w_expect_total;
		result[i].nrOfFurtherModules = internal[i].nrOfFurtherModules;
		result[i].n_a		= internal[i].n_a;
		result[i].n_b		= internal[i].n_b;
		result[i].n		= internal[i].n;
		result[i].label		= internal[i].label;
		result[i].index		= internal[i].index;

		if(internal[i].L->type == DENDRO) {
			result[i].L			= &result[internal[i].L->index];
			result[internal[i].L->index].M	= &result[i];
		}
		else {
			result[i].L			= &leaves[internal[i].L->index];
			leaves[internal[i].L->index].M	= &result[i];
		}

		if(internal[i].R->type == DENDRO) {
			result[i].R			= &result[internal[i].R->index];
			result[internal[i].R->index].M	= &result[i];
		}
		else {
			result[i].R			= &leaves[internal[i].R->index];
			leaves[internal[i].R->index].M	= &result[i];
		}
	}

	return result;
}

// ********************************************************************************************************

// This method merges the two lists leftVertices and rightVertices into one with its elements' x in increasing order.
list* dendro::mergeLists(list* leftVertices, list* rightVertices) {

	list* current;
	list* head;

	if(leftVertices->x < rightVertices->x) {
		head		= leftVertices;
		current		= leftVertices;
		leftVertices	= leftVertices->next;
	}
	else if(leftVertices->x > rightVertices->x) {
		head		= rightVertices;
		current		= rightVertices;
		rightVertices	= rightVertices->next;
	}
	else cout << "!! Error: same vertex" << leftVertices->x << " in both subtrees of internal vertex" << endl;

	while(leftVertices != NULL && rightVertices != NULL) {
		if(leftVertices->x < rightVertices->x) {
			current->next	= leftVertices;
			current		= leftVertices;
			leftVertices	= leftVertices->next;
		}
		else if(leftVertices->x > rightVertices->x) {
			current->next	= rightVertices;
			current		= rightVertices;
			rightVertices	= rightVertices->next;
		}
		else cout << "!! Error: same vertex" << leftVertices->x << " in both subtrees of internal vertex" << endl;
	}

	if(leftVertices != NULL) {
		current->next	= leftVertices;
		current		= NULL;
		leftVertices	= NULL;
	}
	else  {
		current->next	= rightVertices;
		current		= NULL;
		rightVertices	= NULL;
	}

	return head;
}

// ********************************************************************************************************

void dendro::deleteList(list* head){
	if(head != NULL) {
		list *toDelete;
		while(head->next != NULL) {
			toDelete = head;
			head = head->next;
			delete toDelete;
			toDelete = NULL;
		}
		delete head;
		head = NULL;
	}
	return;
}

// ********************************************************************************************************

// Following function writes the order files and module file
// O(n*log(n)) in best case, O(n^2) in worst case (already summed up for all recursive invocations)
list* dendro::recordOrderAndModules(rbtree& reverseNamesLUT, ofstream& orderAFOut, ofstream& orderBFOut, ofstream& modulesFOut, ofstream& nestedModulesFOut, elementd* vertex, const int nrFurtherCompPrevVertex, const int depth) {

	if(vertex->type == DENDRO) {

		list* leftVertices = recordOrderAndModules(reverseNamesLUT, orderAFOut, orderBFOut, modulesFOut, nestedModulesFOut, vertex->L, vertex->nrOfFurtherModules, depth+1);
		list* rightVertices = recordOrderAndModules(reverseNamesLUT, orderAFOut, orderBFOut, modulesFOut, nestedModulesFOut, vertex->R, vertex->nrOfFurtherModules, depth+1);

		list* result = mergeLists(leftVertices, rightVertices);
		list* out = result;

		if(vertex->nrOfFurtherModules <= 0) {

			int i = 1;

			if(!strcmp(task, "computeModulesExperimental")) nestedModulesFOut << depth << "\t";
			if(vertex->nrOfFurtherModules == 0) modulesFOut << depth << "\t";						// write nesting depth of module to module file

			while(out != NULL) {
				while(i < out->x) {
					if(!strcmp(task, "computeModulesExperimental")) nestedModulesFOut << "0\t";
					if(vertex->nrOfFurtherModules == 0) modulesFOut << "0\t";
					i++;
				}
				if(!strcmp(task, "computeModulesExperimental")) nestedModulesFOut << out->x;
				if(vertex->nrOfFurtherModules == 0) modulesFOut << out->x;
				if(out->x != n) {
					if(!strcmp(task, "computeModulesExperimental")) nestedModulesFOut << "\t";
					if(vertex->nrOfFurtherModules == 0) modulesFOut << "\t";
				}
				out = out->next;
				i++;
			}

			while(i <= n) {
				if(!strcmp(task, "computeModulesExperimental")) nestedModulesFOut << "0";
				if(vertex->nrOfFurtherModules == 0) modulesFOut << "0";
				if(i != n) {
					if(!strcmp(task, "computeModulesExperimental")) nestedModulesFOut << "\t";
					if(vertex->nrOfFurtherModules == 0) modulesFOut << "\t";
				}
				i++;
			}

			if(!strcmp(task, "computeModulesExperimental")) nestedModulesFOut << endl;
			if(vertex->nrOfFurtherModules == 0) modulesFOut << endl;
		}

		return result;
	}
	else {
		int thisVerticesLabel = reverseNamesLUT.findItem(vertex->index)->value;

		list* thisVertex	= new list;
		thisVertex->x	= thisVerticesLabel;
		thisVertex->next	= NULL;

		if(vertex->partition == PARTITION_A) {
			orderAFOut << thisVerticesLabel << endl;					// write real(!) index to "order file"
		}
		else if(vertex->partition == PARTITION_B) {
			orderBFOut << thisVerticesLabel << endl;					// write real(!) index to "order file"
		}

		return thisVertex;
	}
}

// ********************************************************************************************************

list* dendro::getInternalVertexIndicesWithinModules() {

	list* head = new list;

	int count = 0;

	for(int i = 0; i < n-1; i++) {
		if(internal[i].nrOfFurtherModules <= 0 && (internal[i].L->type == DENDRO || internal[i].L->type == DENDRO)) {
			head->x		= i;
			count++;
			list* newHead	= new list;
			newHead->next	= head;
			head		= newHead;
		}
	}

	head->x = count;
	return head;
}

// ********************************************************************************************************
// *** public methods *************************************************************************************
// ********************************************************************************************************

void dendro::buildDendrogram(const bool flag_bipartite) {
	if(g == NULL) { cout << "!! Error: cannot build dendrogram without a graph structure.\n"; return; }

	isBipartite = flag_bipartite;

/* the initialization of the dendrogram structure goes like this:
 * 1) we allocate space for the n-1 internal vertices of the dendrogram, and then the n leaf vertices
 * 2) we build a random binary tree structure out of the internal vertices by assigning each
 *    a uniformly random value over [0,1] and then inserting it into the tree according to the 
 *    binary-search rule.
 * 3) next, we make a random permutation of the n leaf vertices and add them to the dendrogram D by
 *    replacing the empty spots in-order
 * 4) then, we compute the path from the root to each leaf and store that in each leaf (this is
 *    prep work for the next step)
 * 5) finally, we compute the values for nL, nR, e (and thus p) and the label for each internal 
 *    vertex by allocating each of the m edges in g to the appropriate internal vertex
 */
	
	// --- Initialization and memory allocation for data structures
	// After allocating the memory for D and G, we need to mark the vertices for G as being
	// non-internal vertices, and then insert them into a random binary tree structure.
	// For simplicity, we make the first internal vertex in the array the root.
	
	bool flag_debug	= true;
	n_a		= g->getNumAVertices();		// number of vertices in Partition A
	n_b		= g->getNumBVertices();		// number of vertices in Partition B
	n		= g->getNumVertices();		// total number of vertices in graph
	sumEdgeWeight	= g->getSumEdgeWeight();	// total sum of edge weights

	leaf		= new elementd [n];		// allocate memory for G, O(n)
	internal	= new elementd [n-1];		// allocate memory for D, O(n)
	d		= new interns(n-2);      	// allocate memory for internal edges of D, O(n)

	for (int i=0; i<n; i++) {			// initialize leaf vertices
		leaf[i].type = GRAPH;
		if(i < n_a) {
			leaf[i].partition	= PARTITION_A;
			leaf[i].n_a		= 1;
			leaf[i].n_b		= 0;
		}
		else {
			leaf[i].partition	= PARTITION_B;
			leaf[i].n_a		= 0;
			leaf[i].n_b		= 1;
		}

		leaf[i].e_w_total	= 0;
		leaf[i].label		= i;
		leaf[i].index		= i;
		leaf[i].n		= 1;
		leaf[i].nrOfFurtherModules = -1;
	}

	if(flag_debug) { cout << ">> dendro: allocated memory for internal and leaf arrays" << endl; }
	root		= &internal[0];		// initialize internal vertices
	root->label	= 0;
	root->index	= 0;
	root->p		= mtr.randExc();
	for (int i=1; i<(n-1); i++) {		// insert remaining internal vertices, O(n log n)
		internal[i].label	= i;
		internal[i].index	= i;
		internal[i].p		= mtr.randExc();
		binarySearchInsert(root, &internal[i]);
	}
	if(flag_debug) { cout << ">> dendro: inserted internal vertices into random binary tree" << endl; }

	// --- Hang leaf vertices off end of dendrogram O(n log n)
	// To impose this random hierarchical relationship on G, we first take a random permutation
	// of the leaf vertices and then replace the NULLs at the bottom of the tree in-order with
	// the leafs. As a hack to ensure that we can find the leafs later using a binary search,
	// we assign each of them the p value of their parent, perturbed slightly so as to preserve
	// the binary search property.

	block* array; array = new block [n];
	for (int i=0; i<n; i++) { array[i].x = mtr.randExc();  array[i].y = i; }
	QsortMain(array, 0, n-1);

	int k=0;						// replace NULLs with leaf vertices, and
	for (int i=0; i<(n-1); i++) {				// maintain binary search property, O(n)
		if(internal[i].L == NULL) {
			internal[i].L = &leaf[array[k].y];
			leaf[array[k].y].M = &internal[i];
			leaf[array[k++].y].p = internal[i].p - 0.0000000000001;
		}
		if(internal[i].R == NULL) {
			internal[i].R = &leaf[array[k].y];
			leaf[array[k].y].M = &internal[i];
			leaf[array[k++].y].p = internal[i].p + 0.0000000000001;
		}
	}
	delete [] array;
	array = NULL;

	if(flag_debug) { cout << ">> dendro: replaced NULLs in bin-tree with leaf vertices" << endl; }

	// --- Compute the path from root -> leaf for each leaf O(n log n)
	// Using the binary search property, we can find each leaf vertex in O(log n) time. The
	// binarySearchFind() method returns the list of internal vertex indices that the search
	// crossed, in the order of root -> ... -> leaf, for use in the subsequent few operations.
	
	if(paths != NULL) {
		list *curr, *prev;

		for (int i = 0; i < n; i++) {
			curr = paths[i];
			while (curr != NULL) {
				prev = curr;
				curr = curr->next;
				delete prev;
				prev = NULL;
			}
			paths[i] = NULL;
		}

		delete [] paths;
	}
	
	paths = NULL;
	paths = new list* [n];
	for (int i=0; i<n; i++) { paths[i] = binarySearchFind(leaf[i].p); }

	if(flag_debug) { cout << ">> dendro: computed paths from root to leafs" << endl; }

	// --- Count e for each internal vertex O(m)
	// To count the number of edges that span the L and R subtrees for each internal vertex and
	// the sum of their weights e_w we use the path information we just computed.
	// Then, we loop over all edges in G and find the common ancestor in D of the two endpoints
	// and increment that internal vertex's e, e_w and e_w_expect.
	// This process takes O(m) time because in a roughly balanced binary
	// tree (given by our random dendrogram), the vast majority of vertices take basically
	// constant time to find their common ancestor.
	//
	// Furthermore, to compute e_w_expect, i.e. the sum of the expected weights of all possible edges
	// between the L and R subtrees of each internal vertex we have to do a computation analogous
	// to the previous one but now for all(!) possible edges.
	//
	// Note that because our adjacency list is symmetric, we overcount each e, e_w and e_w_expect
	// by a factor of 2, so we need to correct this after.

	elementd* ancestor; edge* curr;

	for (int i=0; i<(n-1); i++) {
		internal[i].e = 0;
		internal[i].e_w = 0;
		internal[i].e_w_total = 0;
		internal[i].e_w_expect = 0.0;
		internal[i].label = -1;
		internal[i].nrOfFurtherModules = -1;
	}

	for (int i=0; i<n; i++) {

		curr = g->getNeighborList(i);
		while (curr != NULL) {
			ancestor		= findCommonAncestor(paths, i, curr->x);
			ancestor->e		+= 1;
			ancestor->e_w		+= curr->weight;
			curr			= curr->next;
		}

		// for computing e_w_expect, i.e. the sum of expected weights of the possible edges between each internal
		// vertex's subtrees L and R, we have to iterate over all(!) possible edges
		if(!strcmp(task, "computeModules")) {
			for(int j = 0; j < n; j++) {
				if(i != j) {
					ancestor		= findCommonAncestor(paths, i, j);
					ancestor->e_w_expect	+= g->getExpectedEdgeWeight(i, j);
				}
			}
		}
	}

	for (int i=0; i<(n-1); i++) {
		internal[i].e		/= 2;
		internal[i].e_w		/= 2;
		internal[i].e_w_expect	/= 2.0;
	}

	if(flag_debug) { cout << ">> dendro: finished common ancestor computation" << endl; }

	// --- Count n_a, n_b for each internal vertex O(n log n)
	// To tabulate the number of leafs of each partition in each subtree rooted at an internal vertex,
	// we use the path information computed above.

	bool currentleafIsInPartitionA;
	for (int i=0; i<n; i++) {
		ancestor = &leaf[i];

		if(ancestor->partition == PARTITION_A) {
			currentleafIsInPartitionA = true;
		}
		else {
			currentleafIsInPartitionA = false;
		}
		
		ancestor = ancestor->M;
		while (ancestor != NULL) {
			if(currentleafIsInPartitionA) {
				ancestor->n_a++;
			}
			else {
				ancestor->n_b++;
			}
			ancestor->n++;
			ancestor = ancestor->M;
		}

	}
	
	if(flag_debug) { cout << ">> dendro: computed subtree sizes" << endl; }
	
	// --- Label all internal vertices O(n log n)
	// We want to label each internal vertex with the smallest leaf index of its children.
	// This will allow us to collapse many leaf-orderings into a single dendrogram structure
	// that is independent of child-exhanges (since these have no impact on the likelihood
	// of the hierarchical structure). To do this, we loop over the leaf vertices from 
	// smallest to largest and walk along that leaf's path from the root. If we find an 
	// unlabeled internal vertex, then we mark it with this leaf's index.

	for (int i=0; i<n; i++) {
		ancestor = &leaf[i];
		while (ancestor != NULL) {
			if(ancestor->label == -1 || ancestor->label > leaf[i].label) { ancestor->label = leaf[i].label; }
			ancestor = ancestor->M;
		}
	}
	
	if(flag_debug) { cout << ">> dendro: labeled all internal vertices" << endl; }
	
	// --- Exchange children to enforce order-property O(n)
	// We state that the order-property requires that an internal vertex's label is the
	// smallest index of its left subtree. The dendrogram so far doesn't reflect this, so we
	// need to step through each internal vertex and make that adjustment (swapping nL and nR
	// if we make a change).
	
	elementd *tempe;
	for (int i=0; i<(n-1); i++) {
		if(internal[i].L->label > internal[i].label) {
			tempe          = internal[i].L;
			internal[i].L  = internal[i].R;
			internal[i].R  = tempe;
		}
	}
	if(flag_debug) { cout << ">> dendro: enforced order-property" << endl; }

	// --- Tabulate internal dendrogram edges O(n^2)
	// For the MCMC moves later on, we'll need to be able to choose, uniformly at random, an
	// internal edge of the dendrogram to manipulate. There are always n-2 of them, and we can
	// find them simply by scanning across the internal vertices and observing which have children
	// that are also internal vertices. Note: very important that the order property be enforced
	// before this step is taken; otherwise, the internal edges wont reflect the actual dendrogram
	// structure.
	
	for (int i=0; i<(n-1); i++) {
		if(internal[i].L->type == DENDRO) { d->addEdge(i, internal[i].L->index, LEFT); }
		if(internal[i].R->type == DENDRO) { d->addEdge(i, internal[i].R->index, RIGHT); }
	}
	
	if(flag_debug) { cout << ">> dendro: tabulated internal dendrogram edges" << endl; }
	
	// --- Clear memory for paths O(n log n)
	// Now that we're finished using the paths, we need to deallocate them manually.
	
	list *current, *previous;
	for (int i=0; i<n; i++) {
		current = paths[i];
		while (current != NULL) { previous = current;   current = current->next;   delete previous;   previous = NULL; }
		paths[i] = NULL;
	}
	delete [] paths;
	paths = NULL;
	if(flag_debug) { cout << ">> dendro: cleared memory for paths" << endl; }
	
	setTotalEdgeWeight(root);
	setTotalExpectedEdgeWeight(root);
	setNrOfFurtherModules(root, 1, false, false);

	// --- Compute p_i for each internal vertex O(n)
	// Each internal vertex's p_i = e_i / [ (n_aL_i*n_bR_i) + (n_bL_i*n_aR_i) ], and now that we have each of those
	// pieces, we may calculate this value for each internal vertex. Given these, we can then
	// calculate the likelihood of the entire dendrogram structure

	L = 0.0;

	int nL_nR, ei, ew;
	double dL, ew_expect;

	for (int i=0; i<(n-1); i++) {
		nL_nR		= (internal[i].L->n_a*internal[i].R->n_b + internal[i].L->n_b*internal[i].R->n_a);
		ei		= internal[i].e;
		ew		= internal[i].e_w;
		ew_expect	= internal[i].e_w_expect;

		if(isBipartite && nL_nR == 0 && ei > 0) {
			cout << "**** WARNING - violation of bipartite structure" << internal[i].L->n_a << ", " << internal[i].R->n_b << ", " << internal[i].L->n_b << ", " << internal[i].R->n_a << ", " << endl;
			return;
		}
		else {
			if(!strcmp(task, "computeModules")) {
				internal[i].p = 0.0;
				dL = ((double)(ew) - ew_expect) / (double)(sumEdgeWeight);
				if(internal[i].nrOfFurtherModules == 1) {
					dL = -dL;
				}
				else if(ew == 0) {
					dL = (double)(-(n-1));
				}
			}
			else if(!strcmp(task, "computeModulesExperimental")) {
				internal[i].p = 0.0;
				dL = ((double)(ew) - ew_expect) / (double)(sumEdgeWeight) / (double)(nL_nR);
				if(ew == 0) {
					dL = (double)(-(n-1));
				}
			}
		}

		internal[i].Lcont = dL;
		L += dL;
	}
	
	if(flag_debug) {
		cout << ">> dendro: computed likelihood" << endl;
		cout << "\n>> likelihood = " << L << "\n" << endl;
	}
	char pauseme;
	for (int i=0; i<(n-1); i++) {
		if(internal[i].label > internal[i].L->label) {
			tempe          = internal[i].L;
			internal[i].L  = internal[i].R;
			internal[i].R  = tempe;
			cout << "#### WARNING - order property violated by internal[" << i << "] (fixed)" << endl;
			cin >> pauseme;
		}
	}

	// --- Dendrogram is now built
	if(flag_debug) { cout << ">> dendro: build dendrogram complete" << endl; }

	return;
}

// ********************************************************************************************************

double dendro::getLikelihood() {
	return L;
}

// ********************************************************************************************************

bool dendro::importDendrogramStructure(const string in_file) {
	string bracketL, bracketR, sL, sR, sLtype, sRtype, sp, se, sn_a, sn_b, sn;
	int sindex, sLindex, sRindex, snume, snumn_a, snumn_b, snumn;
	double sprob;
	bool safeExit   = true;
	bool flag_debug = true;
	
	ifstream fscan(in_file.c_str(), ios::in);
	fscan >> bracketL >> sindex >> bracketR >> sL >> sLindex >> sLtype >> sR >> sRindex >> sRtype >> sp >> sprob >> se >> snume >> sn_a >> snumn_a >> sn_b >> snumn_b >> sn >> snumn;
	n_a = snumn_a;
	n_b = snumn_b;
	n = snumn;
	fscan.close();
	
	leaf		= new elementd [n];		// allocate memory for G, O(n)
	internal	= new elementd [n-1];		// allocate memory for D, O(n)
	d		= new interns(n-2);		// allocate memory for internal edges of D, O(n)
	for (int i=0; i<n; i++) {			// initialize leaf vertices
		leaf[i].type	= GRAPH;

		if(i < n_a) {
			leaf[i].partition = PARTITION_A;
			leaf[i].n_a = 1;
			leaf[i].n_b = 0;
		}
		else {
			leaf[i].partition = PARTITION_B;
			leaf[i].n_a = 0;
			leaf[i].n_b = 1;
		}
		
		leaf[i].label	= i;
		leaf[i].index	= i;
		leaf[i].n	= 1;
	}
	root = &internal[0];				// initialize internal vertices
	root->label = 0;
	for (int i=1; i<(n-1); i++) { internal[i].index = i; internal[i].label = -1; }
	if(flag_debug) { cout << ">> dendro: allocated memory for internal and leaf arrays" << endl; }
	
	// --- Import basic structure from file O(n)
	ifstream fin(in_file.c_str(), ios::in);
	while (fin >> bracketL >> sindex >> bracketR >> sL >> sLindex >> sLtype >> sR >> sRindex >> sRtype >> sp >> sprob >> se >> snume >> sn_a >> snumn_a >> sn_b >> snumn_b >> sn >> snumn) {
		cout << bracketL << " " << sindex << " " << bracketR << " " << sL << " " << sLindex << " " << sLtype << " " << sR << " " << sRindex << " " << sRtype << " " << sp << " " << sprob << " " << se << " " << snume << " " << sn_a << " " << snumn_a << " " << sn_b << " " << snumn_b << " " << sn << " " << snumn << endl;

		if(sLtype == "(D)") {
			internal[sindex].L = &internal[sLindex]; internal[sLindex].M = &internal[sindex];
		}
		else if(sLtype == "(G)") {
			internal[sindex].L = &leaf[sLindex];     leaf[sLindex].M     = &internal[sindex];
		}
		else {
			cout << "!! Error: " << bracketL << sindex << bracketR << sL << sLindex << sLtype << sR << sRindex << sRtype << sp << sprob << se << snume << sn_a << snumn_a << sn_b << snumn_b << sn << snumn << endl;
			safeExit = false;
			break;
		}

		if(sRtype == "(D)") {
			internal[sindex].R = &internal[sRindex];
			internal[sRindex].M = &internal[sindex];
		}
		else if(sRtype == "(G)") {
			internal[sindex].R = &leaf[sRindex];
			leaf[sRindex].M = &internal[sindex];
		} else {
			cout << "!! Error: " << bracketL << sindex << bracketR << sL << sLindex << sLtype << sR << sRindex << sRtype << sp << sprob << se << snume << sn_a << snumn_a << sn_b << snumn_b << sn << snumn << endl;
			safeExit = false;
			break;
		}

		internal[sindex].p     = sprob;
		if(sprob < 0.0 || sprob > 1.0) {
			cout << "!! Error: " << bracketL << sindex << bracketR << sL << sLindex << sLtype << sR << sRindex << sRtype << sp << sprob << se << snume << sn_a << snumn_a << sn_b << snumn_b << sn << snumn << endl;
			safeExit = false;
			break;
		}

		internal[sindex].e     = snume;
		internal[sindex].n_a     = snumn_a;
		internal[sindex].n_b     = snumn_b;
		internal[sindex].n     = snumn;
		internal[sindex].index = sindex;
	}
	fin.close();
	if(!safeExit) { return false; }
	if(flag_debug) { cout << ">> dendro: imported basic structure" << endl; }

	// --- Label all internal vertices O(n log n)
	elementd* curr;
	for (int i=0; i<n; i++) {
		curr = &leaf[i];
		while (curr != NULL) { if(curr->label == -1 || curr->label > leaf[i].label) { curr->label = leaf[i].label; } curr = curr->M; }
	}
	if(flag_debug) { cout << ">> dendro: labeled all internal vertices" << endl; }
	
	// --- Exchange children to enforce order-property O(n)
	elementd *tempe;
	for (int i=0; i<(n-1); i++) {
		if(internal[i].L->label > internal[i].label) {
			tempe          = internal[i].L;
			internal[i].L  = internal[i].R;
			internal[i].R  = tempe;
		}
	}
	if(flag_debug) { cout << ">> dendro: enforced order-property" << endl; }
	
	// --- Tabulate internal dendrogram edges O(n)
	for (int i=0; i<(n-1); i++) {
		if(internal[i].L->type == DENDRO) { d->addEdge(i, internal[i].L->index, LEFT);  }
		if(internal[i].R->type == DENDRO) { d->addEdge(i, internal[i].R->index, RIGHT); }
	}
	if(flag_debug) { cout << ">> dendro: tabulated internal dendrogram edges" << endl; }

	// --- Compute p_i for each internal vertex O(n)
	// Each internal vertex's p_i = e_i / [ (n_aL_i*n_bR_i) + (n_bL_i*n_aR_i) ], and now that we have each of those
	// pieces, we may calculate this value for each internal vertex. Given these, we can then
	// calculate the likelihood of the entire dendrogram structure
	// L = \sum_{i=1}^{n} ( p_i )

	L = 0.0;

	int nL_nR, ei, ew;
	double dL, ew_expect;

	for (int i=0; i<(n-1); i++) {
		nL_nR		= (internal[i].L->n_a*internal[i].R->n_b + internal[i].L->n_b*internal[i].R->n_a);
		ei		= internal[i].e;
		ew		= internal[i].e_w;
		ew_expect	= internal[i].e_w_expect;

		if(isBipartite && nL_nR == 0 && ei > 0) {
			cout << "**** WARNING - violation of bipartite structure" << internal[i].L->n_a << ", " << internal[i].R->n_b << ", " << internal[i].L->n_b << ", " << internal[i].R->n_a << ", " << endl;
			return false;
		}
		else {
			if(!strcmp(task, "computeModules")) {
				internal[i].p = 0.0;
				dL = ((double)(ew) - ew_expect) / (double)(sumEdgeWeight);
				if(internal[i].nrOfFurtherModules == 1) {
					dL = -dL;
				}
				else if(ew == 0) {
					dL = (double)(-(n-1));
				}
			}
			else if(!strcmp(task, "computeModulesExperimental")) {
				internal[i].p = 0.0;
				dL = ((double)(ew) - ew_expect) / (double)(sumEdgeWeight) / (double)(nL_nR);
				if(ew == 0) {
					dL = (double)(-(n-1));
				}
			}
		}

		internal[i].Lcont = dL;
		L += dL;
	}

	if(flag_debug) {
		cout << ">> dendro: computed likelihood" << endl;
		cout << "\n>> likelihood = " << L << "\n" << endl;
	}
	
	// --- Dendrogram is now built
	if(flag_debug) { cout << ">> dendro: build dendrogram complete" << endl; }
	
	return true;
}

// ********************************************************************************************************

dendro* dendro::deepCopy() {
	dendro* bestDendro = new dendro(task);
	if(!bestDendro->setValues(n_a, n_b, sumEdgeWeight, L, d, this)) {
		cout << "!! Error: failed to copy current dendrogram to best dendrogram" << endl;
		delete bestDendro;
		return NULL;
	}
	return bestDendro;
}

// ********************************************************************************************************

bool dendro::monteCarloMove(double& delta, bool& ftaken, const double T, const double bestL) {

	// A single MC move begins with the selection of a random internal edge (a,b) of the
	// dendrogram. This also determines the three subtrees i, j, k that we will rearrange,
	// and we choose uniformly from among the options.
	// 
	// if(a,b) is a left-edge, then we have ((i,j),k), and moves
	// ((i,j),k) -> ((i,k),j)								(alpha move)
	//           -> (i,(j,k)) + enforce order-property for (j,k)				(beta move)
	// 
	// if(a,b) is a right-edge, then we have (i,(j,k)), and moves
	// (i,(j,k)) -> ((i,k),j)								(alpha move)
	//           -> ((i,j),k)								(beta move)
	// 
	// For each of these moves, we need to know what the change in likelihood will be, so
	// that we can determine with what probability we execute the move.
	
	bool		flag_debug = true;
	elementd	*temp, *tempe;
	ipair		*tempPair;
	edgeCountTriple*	ect;
	int		x, y, e_x, e_y, e_w_x, e_w_y, n_a_i, n_b_i, n_a_j, n_b_j, n_a_k, n_b_k, n_x, n_y;
	short int	t;
	double		p_x, p_y, Lcont_x, Lcont_y, dLcont, dF, e_w_expect_x, e_w_expect_y, dLcontInvolvedVertices;
	
	// The remainder of the code executes a single MCMC move, where we sample the dendrograms 
	// proportionally to their likelihoods (i.e., temperature=1, if you're comparing it to the
	// usual MCMC framework). 
	delta    = 0.0;
	ftaken   = false;

	tempPair = d->getRandomEdge();		// returns address; do not deallocate

	x        = tempPair->x;			// copy contents of referenced random edge
	y        = tempPair->y;			// into local variables
	t        = tempPair->t;
	
	if(flag_debug) {
		if(t == LEFT || t == RIGHT) {}
		else {
			cout << " bad edge (i,j,k)" << endl;
		}
	}
	
	if(t == LEFT) {
		if(mtr.randExc() < 0.5) {									// LEFT ALPHA move: ((i,j),k) -> ((i,k),j)

			// We need to calculate the change in the likelihood (dLcont) that would result from
			// this move. Most of the information needed to do this is already available,
			// the exception being e_ik, the number of edges that span the i and k subtrees.
			// I use a slow algorithm O(n) to do this, since I don't know of a better way at
			// this point. (After several attempts to find a faster method, no luck.)

			n_a_i = internal[y].L->n_a;
			n_b_i = internal[y].L->n_b;
			n_a_j = internal[y].R->n_a;
			n_b_j = internal[y].R->n_b;
			n_a_k = internal[x].R->n_a;
			n_b_k = internal[x].R->n_b;

			// Recalculation of values of lower subtree's root y
			n_y		= n_a_i * n_b_k + n_b_i * n_a_k;
			ect		= computeEdgeCount(internal[y].L->index, internal[y].L->type, internal[x].R->index, internal[x].R->type);
			e_y		= ect->e;												// e_ik
			e_w_y		= ect->e_w;
			e_w_expect_y	= ect->e_w_expect;

			delete ect;
			ect		= NULL;

			if(isBipartite && n_y == 0 && e_y > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1) {
						Lcont_y = -Lcont_y;
					}
					else if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight) / (double)(n_y);
					if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
			}

			// Recalculation of values of upper subtree's root x
			n_x		= (n_a_i+n_a_k)*n_b_j + (n_b_i+n_b_k)*n_a_j;
			e_x		= internal[x].e + internal[y].e - e_y;									// e_yj
			e_w_x		= internal[x].e_w + internal[y].e_w - e_w_y;
			e_w_expect_x	= internal[x].e_w_expect + internal[y].e_w_expect - e_w_expect_y;

			if(isBipartite && n_x == 0 && e_x > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1 || (internal[x].nrOfFurtherModules + internal[y].nrOfFurtherModules > -2 && internal[y].R->type == DENDRO)) {
						Lcont_x = -Lcont_x;
					}
					else if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight) / (double)(n_x);
					if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
			}

			dLcontInvolvedVertices = 0.0;

			if(!strcmp(task, "computeModules")) {

				if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].R->type == DENDRO){
					dLcontInvolvedVertices = computeLcont(internal[y].R);
				}
				else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {
					dLcontInvolvedVertices = computeLcont(internal[x].R);
					if(internal[y].R->type == DENDRO) {
						dLcontInvolvedVertices += computeLcont(internal[y].R);
					}
				}
			}

			dLcont	= (Lcont_x - internal[x].Lcont) + (Lcont_y - internal[y].Lcont) - dLcontInvolvedVertices;
			dF	= (L + dLcont) - bestL;
			if(dLcont > 0.0 || mtr.randExc() < exp(dLcont/T)) {							// make LEFT ALPHA move

				ftaken = true;

				if(!strcmp(task, "computeModules")) {

					if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].R->type == DENDRO) {
						internal[x].nrOfFurtherModules = 1;
						internal[y].nrOfFurtherModules = 0;
						setNrOfFurtherModules(internal[y].R, 1, true, false);
					}
					else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {

						setNrOfFurtherModules(internal[x].R, -1, true, false);

						if(internal[y].R->type == GRAPH) {
							internal[x].nrOfFurtherModules = 0;
							internal[y].nrOfFurtherModules = -1;
						}
						else {
							setNrOfFurtherModules(internal[y].R, 1, true, false);
						}
					}
				}

				d->swapEdges(x, internal[x].R->index, RIGHT, y, internal[y].R->index, RIGHT);

				temp			= internal[x].R;					// - swap j and k
				internal[x].R		= internal[y].R;					// 
				internal[y].R		= temp;							// 
				internal[x].R->M	= &internal[x];						// - adjust parent pointers
				internal[y].R->M	= &internal[y];									// 
				internal[y].n_a		= n_a_i + n_a_k;					// - update n_a for [y]
				internal[y].n_b		= n_b_i + n_b_k;					// - update n_b for [y]
				internal[y].n		= n_a_i + n_a_k + n_b_i + n_b_k;			// - update n for [y]
				internal[x].e		= e_x;							// - update e_i for [x] and [y]
				internal[y].e		= e_y;							// 
				internal[x].e_w		= e_w_x;						// - update e_w_i for [x] and [y]
				internal[y].e_w		= e_w_y;						//
				internal[x].e_w_expect	= e_w_expect_x;						// - update e_w_expect_i for [x] and [y]
				internal[y].e_w_expect	= e_w_expect_y;						//
				internal[x].p		= p_x;							// - update p_i for [x] and [y]
				internal[y].p		= p_y;							// 
				internal[x].Lcont	= Lcont_x;						// - update L_i for [x] and [y]
				internal[y].Lcont	= Lcont_y;						//
														// - order-property maintained
				L			+= dLcont;						// - update Lcont
				delta			= dLcont;						//

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = " << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = " << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = " << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = " << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = " << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = " << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
			}
		}
		else {											// LEFT BETA move:  ((i,j),k) -> (i,(j,k))

			n_a_i = internal[y].L->n_a;
			n_b_i = internal[y].L->n_b;
			n_a_j = internal[y].R->n_a;
			n_b_j = internal[y].R->n_b;
			n_a_k = internal[x].R->n_a;
			n_b_k = internal[x].R->n_b;
			
			// Recalculation of values of lower subtree's root y
			n_y		= n_a_j * n_b_k + n_b_j * n_a_k;
			ect		= computeEdgeCount(internal[y].R->index, internal[y].R->type, internal[x].R->index, internal[x].R->type);
			e_y		= ect->e;												// e_jk
			e_w_y		= ect->e_w;
			e_w_expect_y	= ect->e_w_expect;

			delete ect;
			ect		= NULL;

			if(isBipartite && n_y == 0 && e_y > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1) {
						Lcont_y = -Lcont_y;
					}
					else if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight) / (double)(n_y);
					if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
			}
			
			// Recalculation of values of upper subtree's root x
			n_x		= (n_a_j+n_a_k)*n_b_i + (n_b_j+n_b_k)*n_a_i;
			e_x		= internal[x].e + internal[y].e - e_y;									// e_yj
			e_w_x		= internal[x].e_w + internal[y].e_w - e_w_y;
			e_w_expect_x	= internal[x].e_w_expect + internal[y].e_w_expect - e_w_expect_y;

			if(isBipartite && n_x == 0 && e_x > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1 || (internal[x].nrOfFurtherModules + internal[y].nrOfFurtherModules > -2 && internal[y].L->type == DENDRO)) {
						Lcont_x = -Lcont_x;
					}
					else if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight) / (double)(n_x);
					if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
			}

			dLcontInvolvedVertices = 0.0;

			if(!strcmp(task, "computeModules")) {

				if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].L->type == DENDRO){
					dLcontInvolvedVertices = computeLcont(internal[y].L);
				}
				else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {
					dLcontInvolvedVertices = computeLcont(internal[x].R);
					if(internal[y].L->type == DENDRO) {
						dLcontInvolvedVertices += computeLcont(internal[y].L);
					}
				}
			}

			dLcont = (Lcont_x - internal[x].Lcont) + (Lcont_y - internal[y].Lcont) - dLcontInvolvedVertices;
			dF	= (L + dLcont) - bestL;
			if(dLcont > 0.0 || mtr.randExc() < exp(dLcont/T)) {								// make LEFT BETA move

				ftaken = true;

				if(!strcmp(task, "computeModules")) {

					if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].L->type == DENDRO) {
						internal[x].nrOfFurtherModules = 1;
						internal[y].nrOfFurtherModules = 0;
						setNrOfFurtherModules(internal[y].L, 1, true, false);
					}
					else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {

						setNrOfFurtherModules(internal[x].R, -1, true, false);

						if(internal[y].L->type == GRAPH) {
							internal[x].nrOfFurtherModules = 0;
							internal[y].nrOfFurtherModules = -1;
						}
						else {
							setNrOfFurtherModules(internal[y].L, 1, true, false);
						}
					}
				}

				d->swapEdges(y, internal[y].L->index, LEFT, y, internal[y].R->index, RIGHT);

				temp			= internal[y].L;					// - swap L and R of [y]
				internal[y].L		= internal[y].R;					// 
				internal[y].R		= temp;							// 
				d->swapEdges(x, internal[x].R->index, RIGHT, y,internal[y].R->index, RIGHT);
				temp			= internal[x].R;					// - swap i and k
				internal[x].R		= internal[y].R;					// 
				internal[y].R		= temp;							// 
				internal[x].R->M	= &internal[x];						// - adjust parent pointers
				internal[y].R->M	= &internal[y];						// 
				d->swapEdges(x, internal[x].L->index, LEFT, x, internal[x].R->index, RIGHT);
				temp			= internal[x].L;					// - swap L and R of [x]
				internal[x].L		= internal[x].R;					// 
				internal[x].R		= temp;							// 
				internal[y].n_a		= n_a_j + n_a_k;					// - update n_a
				internal[y].n_b		= n_b_j + n_b_k;					// - update n_b
				internal[y].n		= n_a_j + n_a_k + n_b_j + n_b_k;			// - update n
				internal[x].e		= e_x;							// - update e_i
				internal[y].e		= e_y;							// 
				internal[x].e_w		= e_w_x;						// - update e_w_i
				internal[y].e_w		= e_w_y;						//
				internal[x].e_w_expect	= e_w_expect_x;						// - update e_w_expect_i for [x] and [y]
				internal[y].e_w_expect	= e_w_expect_y;						//
				internal[x].p		= p_x;							// - update p_i
				internal[y].p		= p_y;							// 
				internal[x].Lcont	= Lcont_x;							// - update Lcont_i
				internal[y].Lcont	= Lcont_y;							//
				if(internal[y].R->label < internal[y].L->label) {
					d->swapEdges(y, internal[y].L->index, LEFT, y, internal[y].R->index, RIGHT);
					temp		= internal[y].L;    					// - enforce order-property if necessary
					internal[y].L	= internal[y].R;    					// 
					internal[y].R	= temp;							// 
				}										// 
				internal[y].label	= internal[y].L->label;					// 
				L			+= dLcont;						// - update Lcont
				delta			= dLcont;						//

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = " << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = " << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = " << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = " << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = " << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = " << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
			}
		}
	} else {												// right-edge: t == RIGHT
		if(mtr.randExc() < 0.5) {									// RIGHT alpha move: (i,(j,k)) -> ((i,k),j)

			n_a_i = internal[x].L->n_a;
			n_b_i = internal[x].L->n_b;
			n_a_j = internal[y].L->n_a;
			n_b_j = internal[y].L->n_b;
			n_a_k = internal[y].R->n_a;
			n_b_k = internal[y].R->n_b;
			
			// Recalculation of values of lower subtree's root y
			n_y		= n_a_i * n_b_k + n_b_i * n_a_k;
			ect		= computeEdgeCount(internal[x].L->index, internal[x].L->type, internal[y].R->index, internal[y].R->type);
			e_y		= ect->e;												// e_ik
			e_w_y		= ect->e_w;
			e_w_expect_y	= ect->e_w_expect;

			delete ect;
			ect		= NULL;

			if(isBipartite && n_y == 0 && e_y > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1) {
						Lcont_y = -Lcont_y;
					}
					else if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight) / (double)(n_y);
					if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
			}
			
			// Recalculation of values of upper subtree's root x
			n_x		= (n_a_i+n_a_k)*n_b_j + (n_b_i+n_b_k)*n_a_j;
			e_x		= internal[x].e + internal[y].e - e_y;									// e_yj
			e_w_x		= internal[x].e_w + internal[y].e_w - e_w_y;
			e_w_expect_x	= internal[x].e_w_expect + internal[y].e_w_expect - e_w_expect_y;

			if(isBipartite && n_x == 0 && e_x > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1 || (internal[x].nrOfFurtherModules + internal[y].nrOfFurtherModules > -2 && internal[y].L->type == DENDRO)) {
						Lcont_x = -Lcont_x;
					}
					else if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight) / (double)(n_x);
					if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
			}

			dLcontInvolvedVertices = 0.0;


			if(!strcmp(task, "computeModules")) {

				if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].L->type == DENDRO){
					dLcontInvolvedVertices = computeLcont(internal[y].L);
				}
				else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {
					dLcontInvolvedVertices = computeLcont(internal[x].L);
					if(internal[y].L->type == DENDRO) {
						dLcontInvolvedVertices += computeLcont(internal[y].L);
					}
				}
			}

			dLcont = (Lcont_x - internal[x].Lcont) + (Lcont_y - internal[y].Lcont) - dLcontInvolvedVertices;
			dF	= (L + dLcont) - bestL;
			if(dLcont > 0.0 || mtr.randExc() < exp(dLcont/T)) {								// make RIGHT ALPHA move

				ftaken = true;

				if(!strcmp(task, "computeModules")) {

					if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].L->type == DENDRO) {
						internal[x].nrOfFurtherModules = 1;
						internal[y].nrOfFurtherModules = 0;
						setNrOfFurtherModules(internal[y].L, 1, true, false);
					}
					else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {

						setNrOfFurtherModules(internal[x].L, -1, true, false);

						if(internal[y].L->type == GRAPH) {
							internal[x].nrOfFurtherModules = 0;
							internal[y].nrOfFurtherModules = -1;
						}
						else {
							setNrOfFurtherModules(internal[y].L, 1, true, false);
						}
					}
				}

				d->swapEdges(x, internal[x].L->index, LEFT, x, internal[x].R->index, RIGHT);

				temp			= internal[x].L;					// - swap L and R of [x]
				internal[x].L		= internal[x].R;					// 
				internal[x].R		= temp;							// 
				d->swapEdges(y, internal[y].L->index, LEFT, x, internal[x].R->index, RIGHT);
				temp			= internal[y].L;					// - swap i and j
				internal[y].L		= internal[x].R;					// 
				internal[x].R		= temp;							// 
				internal[x].R->M	= &internal[x];						// - adjust parent pointers
				internal[y].L->M	= &internal[y];						// 
				internal[y].n_a		= n_a_i + n_a_k;					// - update n_a
				internal[y].n_b		= n_b_i + n_b_k;					// - update n_b
				internal[y].n		= n_a_i + n_a_k + n_b_i + n_b_k;			// - update n
				internal[x].e		= e_x;							// - update e_i
				internal[y].e		= e_y;							// 
				internal[x].e_w		= e_w_x;						// - update e_w_i
				internal[y].e_w		= e_w_y;						//
				internal[x].e_w_expect	= e_w_expect_x;						// - update e_w_expect_i for [x] and [y]
				internal[y].e_w_expect	= e_w_expect_y;						//
				internal[x].p		= p_x;							// - update p_i
				internal[y].p		= p_y;							// 
				internal[x].Lcont	= Lcont_x;							// - update Lcont_i
				internal[y].Lcont	= Lcont_y;							//
				internal[y].label	= internal[x].label;					// - update order property
				L			+= dLcont;						// - update Lcont
				delta			= dLcont;						// 

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = " << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = " << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = " << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = " << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = " << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = " << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
			}
		} else {											// RIGHT beta move:  (i,(j,k)) -> ((i,j),k)

			n_a_i = internal[x].L->n_a;
			n_b_i = internal[x].L->n_b;
			n_a_j = internal[y].L->n_a;
			n_b_j = internal[y].L->n_b;
			n_a_k = internal[y].R->n_a;
			n_b_k = internal[y].R->n_b;
			
			// Recalculation of values of lower subtree's root y
			n_y		= n_a_i * n_b_j + n_b_i * n_a_j;
			ect		= computeEdgeCount(internal[x].L->index, internal[x].L->type, internal[y].L->index, internal[y].L->type);
			e_y		= ect->e;													// e_ij
			e_w_y		= ect->e_w;
			e_w_expect_y	= ect->e_w_expect;

			delete ect;
			ect		= NULL;

			if(isBipartite && n_y == 0 && e_y > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1) {
						Lcont_y = -Lcont_y;
					}
					else if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_y = ((double)(e_w_y) - (e_w_expect_y)) / (double)(sumEdgeWeight) / (double)(n_y);
					if(e_w_y == 0) {
						Lcont_y = (double)(-(n-1));
					}
				}
			}
			
			// Recalculation of values of upper subtree's root x
			n_x		= (n_a_i+n_a_j)*n_b_k + (n_b_i+n_b_j)*n_a_k;
			e_x		= internal[x].e + internal[y].e - e_y;										// e_yk
			e_w_x		= internal[x].e_w + internal[y].e_w - e_w_y;
			e_w_expect_x	= internal[x].e_w_expect + internal[y].e_w_expect - e_w_expect_y;

			if(isBipartite && n_x == 0 && e_x > 0) {
				cout << "**** WARNING - violation of bipartite structure" << endl;
				return false;
			}
			else {
				if(!strcmp(task, "computeModules")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight);
					if(internal[y].nrOfFurtherModules == 1 || (internal[x].nrOfFurtherModules + internal[y].nrOfFurtherModules > -2 && internal[y].R->type == DENDRO)) {
						Lcont_x = -Lcont_x;
					}
					else if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
				else if(!strcmp(task, "computeModulesExperimental")) {
					Lcont_x = ((double)(e_w_x) - (e_w_expect_x)) / (double)(sumEdgeWeight) / (double)(n_x);
					if(e_w_x == 0) {
						Lcont_x = (double)(-(n-1));
					}
				}
			}

			dLcontInvolvedVertices = 0.0;

			if(!strcmp(task, "computeModules")) {

				if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].R->type == DENDRO){
					dLcontInvolvedVertices = computeLcont(internal[y].R);
				}
				else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {
					dLcontInvolvedVertices = computeLcont(internal[x].L);
					if(internal[y].R->type == DENDRO) {
						dLcontInvolvedVertices += computeLcont(internal[y].R);
					}
				}
			}

			dLcont = (Lcont_x - internal[x].Lcont) + (Lcont_y - internal[y].Lcont) - dLcontInvolvedVertices;
			dF	= (L + dLcont) - bestL;
			if(dLcont > 0.0 || mtr.randExc() < exp(dLcont/T)) {								// make RIGHT BETA move

				ftaken = true;

				if(!strcmp(task, "computeModules")) {

					if(internal[x].nrOfFurtherModules == 0 && internal[y].nrOfFurtherModules == -1 && internal[y].R->type == DENDRO) {
						internal[x].nrOfFurtherModules = 1;
						internal[y].nrOfFurtherModules = 0;
						setNrOfFurtherModules(internal[y].R, 1, true, false);
					}
					else if(internal[x].nrOfFurtherModules == 1 && internal[y].nrOfFurtherModules == 0) {

						setNrOfFurtherModules(internal[x].L, -1, true, false);

						if(internal[y].R->type == GRAPH) {
							internal[x].nrOfFurtherModules = 0;
							internal[y].nrOfFurtherModules = -1;
						}
						else {
							setNrOfFurtherModules(internal[y].R, 1, true, false);
						}
					}
				}

				d->swapEdges(x, internal[x].L->index, LEFT, x, internal[x].R->index, RIGHT);

				temp			= internal[x].L;					// - swap L and R of [x]
				internal[x].L		= internal[x].R;					// 
				internal[x].R		= temp;							// 
				d->swapEdges(x, internal[x].R->index, RIGHT, y, internal[y].R->index, RIGHT);
				temp			= internal[x].R;					// - swap i and k
				internal[x].R		= internal[y].R;					// 
				internal[y].R		= temp;							// 
				internal[x].R->M	= &internal[x];						// - adjust parent pointers
				internal[y].R->M	= &internal[y];						// 
				d->swapEdges(y, internal[y].L->index, LEFT, y, internal[y].R->index, RIGHT);
				temp			= internal[y].L;					// - swap L and R of [y]
				internal[y].L		= internal[y].R;					// 
				internal[y].R    	= temp;							// 
				internal[y].n_a		= n_a_i + n_a_j;					// - update n_a
				internal[y].n_b		= n_b_i + n_b_j;					// - update n_b
				internal[y].n		= n_a_i + n_a_j + n_b_i + n_b_j;			// - update n
				internal[x].e		= e_x;							// - update e_i
				internal[y].e		= e_y;							// 
				internal[x].e_w		= e_w_x;						// - update e_w_i
				internal[y].e_w		= e_w_y;						// 
				internal[x].e_w_expect	= e_w_expect_x;						// - update e_w_expect_i for [x] and [y]
				internal[y].e_w_expect	= e_w_expect_y;						//
				internal[x].p		= p_x;							// - update p_i
				internal[y].p		= p_y;							//
				internal[x].Lcont	= Lcont_x;							// - update Lcont_i
				internal[y].Lcont	= Lcont_y;							//
				internal[y].label	= internal[x].label;					// - order-property
				L			+= dLcont;						// - update Lcont
				delta			= dLcont;						// 

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label or internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						cout << "**** WARNING - order property violated by internal[" << x << "]" << endl;
						cout << "x    (p = " << internal[x].p << "\te = "   << internal[x].e << "\tnL = " << internal[x].L->n << "\tnR = " << internal[x].R->n << "\tlabel = " << internal[x].label << ")\tinternal[" << x << "]\t(D)" << endl;
						if (internal[x].L->type==GRAPH) { cout << "x->L [" << internal[x].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[x].L->p << "\te = "   << internal[x].L->e << "\tnL = " << internal[x].L->L->n << "\tnR = " << internal[x].L->R->n << "\tlabel = " << internal[x].L->label << ")\tinternal[" << internal[x].L->index << "]\t(D)"<< endl; }
						if (internal[x].R->type==GRAPH) { cout << "x->R [" << internal[x].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[x].R->p << "\te = "   << internal[x].R->e << "\tnL = " << internal[x].R->L->n << "\tnR = " << internal[x].R->R->n << "\tlabel = " << internal[x].R->label << ")\tinternal[" << internal[x].R->index << "]\t(D)"<< endl; }
					}
					if (internal[y].label > internal[y].L->label) {
						cout << "**** WARNING - order property violated by internal[" << y << "]" << endl;
						cout << "y    (p = " << internal[y].p << "\te = "   << internal[y].e << "\tnL = " << internal[y].L->n << "\tnR = " << internal[y].R->n << "\tlabel = " << internal[y].label << ")\tinternal[" << y << "]\t(D)" << endl;
						if (internal[y].L->type==GRAPH) { cout << "y->L [" << internal[y].L->index << "]\t(G)"<< endl; } else { cout << "i->L (p = " << internal[y].L->p << "\te = "   << internal[y].L->e << "\tnL = " << internal[y].L->L->n << "\tnR = " << internal[y].L->R->n << "\tlabel = " << internal[y].L->label << ")\tinternal[" << internal[y].L->index << "]\t(D)"<< endl; }
						if (internal[y].R->type==GRAPH) { cout << "y->R [" << internal[y].R->index << "]\t(G)"<< endl; } else { cout << "i->R (p = " << internal[y].R->p << "\te = "   << internal[y].R->e << "\tnL = " << internal[y].R->L->n << "\tnR = " << internal[y].R->R->n << "\tlabel = " << internal[y].R->label << ")\tinternal[" << internal[y].R->index << "]\t(D)"<< endl; }
					}
					return false;
				}
				if (internal[x].label > internal[x].L->label) { tempe = internal[x].L; internal[x].L  = internal[x].R; internal[x].R  = tempe; }
				if (internal[y].label > internal[y].L->label) { tempe = internal[y].L; internal[y].L  = internal[y].R; internal[y].R  = tempe; }
			}
		}
	}

	if((L >= 0.0 && L-delta < 0.0) || (L < 0.0 && L-delta >= 0.0)) {
		refreshLikelihood();
	}

	return true;
}

// ********************************************************************************************************

void dendro::refreshLikelihood() {							// recalculates the likelihood of the dendrogram structure

	bool flag_debug = false;

	double L_new = 0.0;

	int nL_nR, ei, ew;
	double dL, ew_expect;
	edgeCountTriple* ect;

	for (int i=0; i<(n-1); i++) {
		nL_nR		= internal[i].L->n_a*internal[i].R->n_b + internal[i].L->n_b*internal[i].R->n_a;
		ei		= internal[i].e;
		ew		= internal[i].e_w;

		ect = computeEdgeCount(internal[i].L->index, internal[i].L->type, internal[i].R->index, internal[i].R->type);

		ew_expect	= ect->e_w_expect;

		if(isBipartite && nL_nR == 0 && ei > 0) {
			cout << "**** WARNING - violation of bipartite structure" << internal[i].L->n_a << ", " << internal[i].R->n_b << ", " << internal[i].L->n_b << ", " << internal[i].R->n_a << ", " << endl;
			return;
		}
		else {
			if(!strcmp(task, "computeModules")) {
				internal[i].p = 0.0;
				dL = ((double)(ew) - ew_expect) / (double)(sumEdgeWeight);
				if(internal[i].nrOfFurtherModules == 1) {
					dL = -dL;
				}
				else if(ew == 0) {
					dL = (double)(-(n-1));
				}
			}
			else if(!strcmp(task, "computeModulesExperimental")) {
				internal[i].p = 0.0;
				dL = ((double)(ew) - ew_expect) / (double)(sumEdgeWeight) / (double)(nL_nR);
				if(ew == 0) {
					dL = (double)(-(n-1));
				}
			}
		}

		internal[i].Lcont = dL;
		L_new += dL;
	}

	if(flag_debug) {
		cout << ">> dendro: refreshed likelihoods" << endl;
		cout << ">> floating point mistake: " << L - L_new << endl;
	}

	L = L_new;

	return;
}

// ********************************************************************************************************

void dendro::printDendrogram()      {
	cout << "\nLEAFS = " << n << endl << "# ";
	printSubTree(root);
	return;
}

// ********************************************************************************************************

void dendro::makeRandomGraph() {
	bool flag_debug = true;
	if(flag_debug) {
		cout << ">> dendro: making random graph from dendrogram" << endl;
	}
	if(g != NULL) {
		delete g;
	}
	g = NULL;
//	g = new graph(n_a, n_b);

	list	*curr, *prev; 
	if(paths != NULL) {
		for (int i=0; i<n; i++) {
			curr = paths[i];
			while (curr != NULL) {
				prev = curr;
				curr = curr->next;
				delete prev;
				prev = NULL;
			}
			paths[i] = NULL;
		}
		delete [] paths;
	}
	paths = NULL;
	paths = new list* [n];				// build paths from root O(n d)
	for (int i=0; i<n; i++) {
		paths[i] = reversePathToRoot(i);
	}

	bool aToB = true;

	elementd* commonAncestor;
	for (int i=0; i<n; i++) {			// O((h+d)*n^2) - h: height of D; d: average degree in G
		for (int j=(i+1); j<n; j++) {		// decide neighbors of v_i
			commonAncestor = findCommonAncestor(paths,i,j);
			if(mtr.randExc() < commonAncestor->p and leaf[i].partition != leaf[j].partition) { 
				if(!(g->doesLinkExist(i,j))) { if(!(g->addLink(i,j,1,aToB))) { cout << "!! Error: (" << j << " " << i << ")" << endl; } }
				if(!(g->doesLinkExist(j,i))) { if(!(g->addLink(j,i,1,!aToB))) { cout << "!! Error: (" << j << " " << i << ")" << endl; } }
			}
		}
	}
	g->printPairs();
	
	for (int i=0; i<n; i++) {
		curr = paths[i];
		while (curr != NULL) {
			prev = curr;
			curr = curr->next;
			delete prev;
			prev = NULL;
		}
		paths[i] = NULL;
	}
	delete [] paths;					// delete paths data structure O(n log n)
	paths = NULL;

	return;
}

// ********************************************************************************************************

void dendro::recordGraphStructure(const string out_file) {
	edge* curr;
	string thisName;
	bool flag_debug = true;
	if(flag_debug) { cout << ">> dendro: writing random graph to file" << endl; }
	
	ofstream fout(out_file.c_str(), ios::trunc);
	for (int i=0; i<n; i++) {
		curr     = g->getNeighborList(i);
		thisName = g->getName(i);
		while (curr != NULL) {
			if(thisName=="") { fout << i << "\t" << curr->x << endl; }
			else {              fout << thisName << "\t" << g->getName(curr->x) << endl; }
			curr = curr->next;
		}
	}
	fout.close();
	
	return;
}

// ********************************************************************************************************

bool dendro::recordOrderAndModules(rbtree& reverseNamesLUT, const string out_orderAFile, const string out_orderBFile, const string out_modulesFile, const string out_nestedModulesFile) {

	// O(n)
	if(!strcmp(task, "computeModulesExperimental")) setBackNrOfFurtherModules(root);

	// For each internal vertex we compute whether there are more modules within its subtrees
	// than the one its two subtrees form together and assign this number to it.
	// Furthermore, after the invocation of setNrOfFurtherModules() nrOfModules represents the number of modules
	// at the lowest level of the binary tree, i.e. those without further nested modules
	// O(n)
	nrOfModules = 0;
	if(!strcmp(task, "computeModulesExperimental")) setNrOfFurtherModules(root);

	// write the vertices' real(!) indices into order files and module file
	ofstream orderAFout(out_orderAFile.c_str(), ios::trunc);
	ofstream orderBFout(out_orderBFile.c_str(), ios::trunc);
	ofstream modulesFOut(out_modulesFile.c_str(), ios::trunc);
	ofstream nestedModulesFOut(out_nestedModulesFile.c_str(), ios::trunc);

	orderAFout << "orderA\n";
	orderBFout << "orderB\n";

	modulesFOut << "depth\t";
	for(int i = 1; i <= n; i++) {
		modulesFOut << "vertex " << i;
		if(i != n) modulesFOut << "\t";
	}
	modulesFOut << endl;

	if(!strcmp(task, "computeModulesExperimental")) {
		nestedModulesFOut << "depth\t";
		for(int i = 1; i <= n; i++) {
			nestedModulesFOut << "vertex " << i;
			if(i != n) nestedModulesFOut << "\t";
		}
		nestedModulesFOut << endl;
	}
	
	// O(n*log(n)) in best case, O(n^2) in worst case
	list* head = recordOrderAndModules(reverseNamesLUT, orderAFout, orderBFout, modulesFOut, nestedModulesFOut, root, root->nrOfFurtherModules+1, 0);
	deleteList(head);

	orderAFout.close();
	orderBFout.close();
	modulesFOut.close();
	nestedModulesFOut.close();
	
	cout << ">> recorded orderA file: " << out_orderAFile << endl;
	cout << ">> recorded orderB file: " << out_orderBFile << endl;
	cout << ">> recorded modules file: " << out_modulesFile << endl;

	return true;
}

// ********************************************************************************************************

void dendro::recordDendrogramStructure(const string out_file) {

	ofstream fout(out_file.c_str(), ios::trunc);
	for (int i=0; i<(n-1); i++) {
		fout << "[ " << i << " ] ";
		fout << "L= " << internal[i].L->index << "\t"; if(internal[i].L->type == DENDRO) { fout << "(D)\t"; } else { fout << "(G)\t"; }
		fout << "R= " << internal[i].R->index << "\t"; if(internal[i].R->type == DENDRO) { fout << "(D)\t"; } else { fout << "(G)\t"; }
		fout << "p= " << internal[i].p << "\t\t";
		fout << "dL= " << internal[i].Lcont << "\t\t";
		fout << "e= " << internal[i].e << "\t";
//		fout << "e_w= " << internal[i].e_w << "\t";
		fout << "n_a= " << internal[i].n_a << "\t";
		fout << "n_b= " << internal[i].n_b << "\t";
//		fout << "expect= " << internal[i].e_w_expect << "\t";
//		fout << "nrComp= " << internal[i].nrOfFurtherModules << "\t";
		fout << "n= " << internal[i].n << endl;
	}

	fout.close();
	
	cout << ">> recorded structure of best dendrogram to file: " << out_file << endl;
	
	return;
}

// ********************************************************************************************************

void dendro::setTask(char* new_task) {
	task = new_task;
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************

#endif
