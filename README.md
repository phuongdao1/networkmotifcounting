# NETWORK MOTIF COUNTING

## COMPILATION:

g++ NetworkMotifCount.cpp -o NetworkMotifCount

change to long double if the number of trees is large but notice that it only works on 64 bit machine.

## INPUT FILES
### Network File Format 

We include a sample network file ExampleNetwork. The first line contains the number of vertices and edges following by the edges together with their confidence score (should be between 0 and 1). All vertex ids start from 0.

### Tree Topology File Format 

Here we enclose all undirected unlabelled tree topologies with 6,7,8 and 9 vertices.

6trees containing all tree topologies of size 6 
7trees containing all tree topologies of size 7
...

The first number in each file is the number of vertices that all the trees in this file contain. Each line following describe a topology of a tree. The format of the tree topology is described through this following example:

0 1 2 3 1 5

There is an edge between vertex 1 and 2, vertex 2 and 3, vertex 3 and 4, vertex 1 and 5, vertex 5 and 6. This is a simple path 6-5-1-2-3-4

## USAGE

./NetworkMotifCount 6trees 1 ExampleNetwork

First parameter is the file containing all tree topoligies, followed by the id of the tree we would like to count and the file describing the network topology.
