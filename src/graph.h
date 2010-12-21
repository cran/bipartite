// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// graph.h - graph data structure for hierarchical random graphs
// Copyright (C) 2005-2008 Aaron Clauset
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
// Created      : fitHRG:	8 November 2005
//		: consensusHRG:	8 November 2005
//		: predictHRG:	21 June 2006
//
// Modified     : fitHRG:	23 December 2007 (cleaned up for public consumption)
//		: consensusHRG:	23 December 2007 (cleaned up for public consumption)
//		: predictHRG:	23 December 2007 (cleaned up for public consumption)
//
// Modified by Rouven Strauss:
//		March 13, 2010:		merged graph.h files of fitHRG, consensusHRG and predictHRG into one
//					merged graph_simp.h into this one
//		March 18, 2010:		modified functions for usage with quantitative networks
//		March 21, 2010:		added datastructure for margin totals likelihood
//
// ****************************************************************************************************
// 
// Graph data structure for maximum likelihood hrgs. The basic structure is an adjacency list of
// edges; however, many additional pieces of metadata are stored as well. Each vertex stores its
// external name and its degree. Each edge stores a histogram of the probabilities assigned to it 
// by the dendrogram structure. Generally, edges are directional, an adjacency (i,j) should also be 
// stored as (j,i) in order to maintain the undirectedness of the graph.
// 
// ****************************************************************************************************

#if !defined(graph_INCLUDED)
#define graph_INCLUDED

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"
#include "math.h"

#include "rbtree.h"

using namespace std;

// ******** Basic Structures ******************************************************************************

#if !defined(edge_INCLUDED)
#define edge_INCLUDED
class edge {
public:

	edge();
	~edge();
	
	int	x;					// index of edge terminator
	int	weight;					// weight of edge
	edge*	next;					// pointer to next elementd
};

edge::edge()  {
	x		= -1;
	weight		= 1;
	next		= NULL;
}

edge::~edge() {

}

#endif

#if !defined(vert_INCLUDED)
#define vert_INCLUDED
class vert {
public:
	vert();
	~vert();

	string	name;					// (external) name of vertex
	int	degree;					// degree of this vertex

};

vert::vert()  {
	name = "";
	degree = 0;

}

vert::~vert() {}
#endif

// ******** Graph Class with Edge Statistics *************************************************************

class graph {
public:
	graph(const int, const int, const bool);
	~graph();

	bool	addLink(const int, const int, const int, const bool);			// add edge (i,j) with weight to graph
	bool	doesLinkExist(const int, const int);					// true if edge (i,j) is already in graph
	int	getEdgeWeight(const int, const int);					// returns weight of edge (i,j)
	int	getDegree(const int);							// returns degree of vertex i
	string	getName(const int);							// returns name of vertex i
	edge*	getNeighborList(const int);						// returns edge list of vertex i
	int	getNumLinks();								// returns m
	int	getNumAVertices();								// returns n_a
	int	getNumBVertices();								// returns n_b
	int	getNumVertices();								// returns n
	int	getSumEdgeWeight();							// returns sumEdgeWeight
	double	getExpectedEdgeWeight(const int, const int);				// returns expected value of edge weight between i and j;
	void	computeExpectedEdgeWeight();						// computes expected edge weights for all i, j
	void	printPairs();								// prints all edges in graph
	bool	setName(const int, const string);					// set name of vertex i
	
private:
	vert*		vertices;			// list of vertices
	edge**		vertexLink;		// linked list of neighbors to vertex
	edge**		vertexLinkTail;		// pointers to tail of neighbor list
	int		n_a;			// number of A vertices
	int		n_b;			// number of B vertices
	int		n;			// number of vertices
	int		sumEdgeWeight;		// total sum of edge weights
	int		m;			// number of directed edges
	bool		isBipartite;		// bool indicating whether graph is bipartite

	// needed only for fitHRG
	double**	expectedEdgeWeight;	// array containing expected weights of edges computed by following formula:
						// expectedEdgeWeight[i][j] = sum(edge_weight_i_b) * sum(edge_weight_a_j) / sumEdgeWeight, 
						// where b B vertix, a A vertix
	int*		marginTotal;		// contains the margin totals needed for computing expectedEdgeWeight
};

// ******** Constructor / Destructor **********************************************************************

graph::graph(const int sizeOfA, const int sizeOfB, const bool flag_bipartite)  {
	n_a			= sizeOfA;
	n_b			= sizeOfB;
	n			= sizeOfA + sizeOfB;
	sumEdgeWeight		= 0;
	m			= 0;
	isBipartite		= flag_bipartite;
	vertices			= new vert  [n];
	vertexLink		= new edge* [n];
	vertexLinkTail		= new edge* [n];

	for(int i = 0; i < n; i++) {
		vertexLink[i] = NULL;
	}

	// needed only for fitHRG
	marginTotal	= new int[n];
	for(int i = 0; i < n; i++) {
		marginTotal[i] = 0;
	}

	expectedEdgeWeight	= new double*[n_a];
	for(int i = 0; i < n_a; i++) {
		if(isBipartite) {
			expectedEdgeWeight[i] = new double[n_b];
			for(int j = 0; j < n_b; j++) {
				expectedEdgeWeight[i][j] = 0.0;
			}
		}
		else {
			expectedEdgeWeight[i] = new double[n_a];
			for(int j = 0; j < n_a; j++) {
				expectedEdgeWeight[i][j] = 0.0;
			}
		}
	}
}

graph::~graph() {
	edge *curr, *prev;
	for (int i=0; i<n; i++) {
		curr = vertexLink[i];
		while (curr != NULL) {
			prev = curr;
			curr = curr->next;
			delete prev;						// deletes edge histogram, too
		}
	}

	curr = NULL;
	prev = NULL;

	delete [] vertexLink;	vertexLink	= NULL;
	delete [] vertexLinkTail;	vertexLinkTail	= NULL;
	delete [] vertices;	vertices		= NULL;				// deletes vertex histogram, too
}

// ********************************************************************************************************

bool graph::addLink(const int i, const int j, const int weight, const bool aToB) {	// adds the directed edge (i,j,weight) to the adjacency list for v_i
	if (i >= 0 and i < n and j >= 0 and j < n and (!isBipartite || ((i < n_a and j >= n_a) or (j < n_a and i >= n_a)))) {

		edge* newedge;

		newedge		= new edge;
		newedge->x	= j;
		newedge->weight	= weight;

		if(aToB) {
			if (!isBipartite || (i < n_a and j >= n_a and j < n)) {
				sumEdgeWeight	+= weight;
				marginTotal[i]	+= weight;
				if(i != j) {
					marginTotal[j]	+= weight;
				}
			}
			else return false;
		}

		if (vertexLink[i] == NULL) {				// first neighbor
			vertexLink[i]	= newedge;
			vertexLinkTail[i]	= newedge;
			vertices[i].degree	= 1;
		}
		else {							// subsequent neighbor
			vertexLinkTail[i]->next = newedge;
			vertexLinkTail[i]       = newedge;
			vertices[i].degree++;
		}
		m++;							// increment edge count
		newedge = NULL;

		return true;
	}
	if(isBipartite && ((i < n_a && j < n_a) || (i >= n_a and j >= n_a))) {
		cout << "!! Error: graph is not bipartite!" << endl;
		cout << "          Do not invoke -bipartite if graph is unipartite." << endl;
	}
	return false;
}

// ********************************************************************************************************

bool graph::doesLinkExist(const int i, const int j) {	// determines if the edge (i,j) already exists in the adjacency list of v_i
	if (i >= 0 and i < n and j >= 0 and j < n and (!isBipartite || (i < n_a and j >= n_a) or (j < n_a and i >= n_a))) {

		edge* curr;

		curr = vertexLink[i];
		while (curr != NULL) {
			if (curr->x == j) {
				curr = NULL;
				return true;
			}
			curr = curr->next;
		}
		curr = NULL;
	}
	return false;
}

// ********************************************************************************************************

int graph::getEdgeWeight(const int i, const int j) {	// returns the weight of edge (i,j)
	if (i >= 0 and i < n and j >= 0 and j < n and (!isBipartite || (i < n_a and j >= n_a) or (j < n_a and i >= n_a))) {
		edge* curr;
		int weight;
		curr = vertexLink[i];
		while (curr != NULL) {
			if (curr->x == j) {
				weight = curr->weight;
				curr = NULL;
				return weight;
			}
			curr = curr->next;
		}
		curr = NULL;
	}
	return -1;
}

// ********************************************************************************************************

int graph::getDegree(const int i) {
	if (i >= 0 and i < n) {
		return vertices[i].degree;
	} 
	else {
		return -1;
	} 
}

// ********************************************************************************************************

string graph::getName(const int i) {
	if (i >= 0 and i < n) {
		return vertices[i].name;
	}
	else {
		return "";
	}
}

// ********************************************************************************************************

// NOTE: The following method returns addresses; deallocation of returned object is dangerous
edge* graph::getNeighborList(const int i) {
	if (i >= 0 and i < n) {
		return vertexLink[i];
	}
	else {
		return NULL;
	}
}

// ********************************************************************************************************

int	graph::getNumLinks()		{ return m; }
int	graph::getNumAVertices()		{ return n_a; }
int	graph::getNumBVertices()		{ return n_b; }
int	graph::getNumVertices()		{ return n; }
int	graph::getSumEdgeWeight()	{ return sumEdgeWeight; }

double	graph::getExpectedEdgeWeight(const int i, const int j) {

	if (i >= 0 and i < n and j >= 0 and j < n) {
		if(isBipartite) {
			if (i < n_a and j >= n_a) return expectedEdgeWeight[i][j-n_a];
			else if (j < n_a and i >= n_a) return expectedEdgeWeight[j][i-n_a];
			else return 0.0;
		}
		else {
			return expectedEdgeWeight[i][j];
		}
	}
	else {
		cout << "Error: trying to get expected edge weight between vertices of which at least one does not exist" << endl;
		return 0.0;
	}
}

void	graph::computeExpectedEdgeWeight() {

	for(int i = 0; i < n_a; i++) {
		if(isBipartite) {
			for(int j = 0; j < n_b; j++) {
				expectedEdgeWeight[i][j] = (double)(marginTotal[i]) * (double)(marginTotal[j+n_a]) / (double)(sumEdgeWeight);
			}
		}
		else {
			for(int j = 0; j < n_a; j++) {
				expectedEdgeWeight[i][j] = (double)(marginTotal[i]) * (double)(marginTotal[j]) / (double)(sumEdgeWeight);
			}
		}
	}

	delete [] marginTotal;
	marginTotal = NULL;
}

// ********************************************************************************************************

void graph::printPairs() {
	edge* curr;
	int edgeCount = 0;
	for (int i=0; i<n; i++) {
		cout << "[" << i << "]\t";
		curr = vertexLink[i];
		while (curr != NULL) {
			cout << curr->x << "\t";
			edgeCount++;
			curr = curr->next;
		}
		cout << "\n";
	}
	cout << edgeCount << " edges total.\n";
	return;
}

// ********************************************************************************************************

bool graph::setName(const int i, const string text) {
	if (i >= 0 and i < n) {
		vertices[i].name = text;
		return true;
	}
	else {
		return false;
	}
}

// ********************************************************************************************************
// ********************************************************************************************************

#endif
