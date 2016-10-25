/*
 * Team Members:
 *    Will Hauber
 *    Matthew McClellan
 *    John Mikolay
 *
 * Date: Nov 6, 2016
 *
 * Assignment: 3
 */

#include <iostream>
#include <queue>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std; 

//---------------//
// Graph Structs //
//---------------//
struct graphEdge
{
	int firstNode;
	int secondNode;
};

struct graphMatrix
{
	int numberOfNodes;
	vector< vector<int>* >* matrix;
};


//-------------------------//
// Graph Method Signatures //
//-------------------------//
graphMatrix generateAdjacencyMatrix(int numNodes, vector<graphEdge> edges);

vector<int> getNeighbors(int node, graphMatrix graphAdjacencyMatrix);
vector<int> getDistancesUsingBfs(int node, graphMatrix graphAdjacencyMatrix);
int getPathDistanceUsingBfs(int node1, int node2, graphMatrix graphAdjacencyMatrix);

void printMatrix(graphMatrix matrix);
void deleteMatrix(graphMatrix matrix);


//------//
// Main //
//------//
int main(int argc, char* argv[])
{
	// Parse command line arguments
	string usageString = "Usage: ./a.out <numNodes> <edge1Node1> <edge1Node2> <edge2Node1 <edge2Node2> ... -1";
	if (argc <= 1)
	{
		cout << usageString << endl;
		return 1;
	}
	
	int numberOfNodes = atoi(argv[1]);
	cout << "Number of nodes: " << numberOfNodes << endl;
	
	vector<graphEdge> edges;	
	for (int ii = 2 /* skip to the edges args */; ii < argc; ii += 2)
	{
		// See if we're at the end
		if (atoi(argv[ii]) < 0)
		{
			break;
		}
		
		// Make sure there are at least two more items
		if (ii + 1 == argc)
		{
			cout << usageString << endl;
			return 2;
		}
		
		// Make sure the next item isn't the end
		if (atoi(argv[ii+1]) < 0)
		{
			cout << usageString << endl;
			return 3;
		}
		
		// Otherwise, add the edge
		graphEdge newEdge;
		newEdge.firstNode = atoi(argv[ii]);
		newEdge.secondNode = atoi(argv[ii+1]);
		edges.push_back(newEdge);
	}
	
	for (int ii = 0; ii < edges.size(); ii++)
	{
		cout << "Edge " << ii << ": " << edges[ii].firstNode << " " << edges[ii].secondNode << endl;
	}
	
	graphMatrix adjacencyMatrix = generateAdjacencyMatrix(numberOfNodes, edges);	
	printMatrix(adjacencyMatrix);
	
	vector<int> distances = getDistancesUsingBfs(4, adjacencyMatrix);
	for (int ii = 0; ii < distances.size(); ii++)
	{
		cout << distances[ii] << " ";
	}
	cout << endl;
	
	cout << getPathDistanceUsingBfs(1, 3, adjacencyMatrix) << endl;;
	
	
	deleteMatrix(adjacencyMatrix);
	return 0;
}

//--------------------------//
// Graph Method Definitions //
//--------------------------//
graphMatrix generateAdjacencyMatrix(int numNodes, vector<graphEdge> edges)
{
	graphMatrix adjacencyMatrix;
	adjacencyMatrix.numberOfNodes = numNodes;
	adjacencyMatrix.matrix = new vector< vector<int>* >();
	
	// Build the matrix row for each node
	for (int currentNode = 0; currentNode < numNodes; currentNode++)
	{
		vector<int>* matrixRow = new vector<int>();
		
		// Initialize the row
		for (int ii = 0; ii < numNodes; ii++)
		{
			matrixRow->push_back(0);
		}
		
		// Find adjancent nodes		
		for (int edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++)
		{
			// If the current node participates in this edge, it is adjacent to the other node on the edge
			graphEdge currentEdge = edges[edgeIndex];
			if (currentEdge.firstNode == currentNode)
			{
				(*matrixRow)[currentEdge.secondNode] = 1;
			}
			if (currentEdge.secondNode == currentNode)
			{
				(*matrixRow)[currentEdge.firstNode] = 1;
			}
		}
		
		// Add the row to the matrix
		adjacencyMatrix.matrix->push_back(matrixRow);
	}
	
	//printMatrix(adjacencyMatrix);
	
	return adjacencyMatrix;
}

vector<int> getNeighbors(int node, graphMatrix graphAdjacencyMatrix)
{
	vector<int>* matrixRowForNode = (*graphAdjacencyMatrix.matrix)[node];
	
	vector<int> neighbors;
	for (int ii = 0; ii < matrixRowForNode->size(); ii++)
	{
		if ((*matrixRowForNode)[ii] == 1)
		{
			neighbors.push_back(ii);
		}
	}
	return neighbors;
}

vector<int> getDistancesUsingBfs(int node, graphMatrix graphAdjacencyMatrix)
{
	// Initialize the distances vector
	vector<int> distances;
	for (int ii = 0; ii < graphAdjacencyMatrix.numberOfNodes; ii++)
	{
		// Assume we can't reach each node to start with
		distances.push_back(-1);
	}
	// We're distance zero from ourself
	distances[node] = 0;
	
	// Initialize the open list
	queue<int> openList;
	openList.push(node);
	
	// Start searching
	while (openList.size() > 0)
	{
		// Remove the next node in the queue
		int currentNode = openList.front();
		openList.pop();
		int currentDistance = distances[currentNode];
				
		// Loop throught this node's neighbors
		vector<int> neighbors = getNeighbors(currentNode, graphAdjacencyMatrix);
		for (int ii = 0; ii < neighbors.size(); ii++)
		{
			int neighbor = neighbors[ii];
			
			// Skip this neighbor if we've already visited it
			if (distances[neighbor] != -1)
			{
				continue;
			}
			
			// Set the distance this neighbor, then add it to the open list
			distances[neighbor] = currentDistance + 1;
			openList.push(neighbor);			
		}
	}
	
	return distances;
}

int getPathDistanceUsingBfs(int node1, int node2, graphMatrix graphAdjacencyMatrix)
{
	vector<int> node1Distances = getDistancesUsingBfs(node1, graphAdjacencyMatrix);
	return node1Distances[node2];
}

void printMatrix(graphMatrix matrixToPrint)
{
	for (int ii = 0; ii < matrixToPrint.matrix->size(); ii++)
	{
		vector<int>* currentRow = (*matrixToPrint.matrix)[ii];
		for (int jj = 0; jj < currentRow->size(); jj++)
		{
			cout << (*currentRow)[jj] << " ";
		}
		
		cout << endl;
	}
}

void deleteMatrix(graphMatrix matrixToDelete)
{
	for (int ii = 0; ii < matrixToDelete.matrix->size(); ii++)
	{
		delete (*matrixToDelete.matrix)[ii];
	}
	delete matrixToDelete.matrix;
}

