# Network Motif Counting

# NETWORK FILE (exampleNetwork):

the first line contain # vertices and # edges
following by the edges together with their confidence score (should be between 0 and 1)
vertex ids start from 0

# TREE TOPOLOGY FILE 

6trees containing all tree topologies of size 6 
7trees...

the first number in each file is the size of all the trees in this file
following by all the undirected unlabelled trees of that size

the format of the tree topology as following:

0 1 2 3 1 5

the vertex with id 2 connects to vertex with id 1
the vertex with id 3 connects to vertex with id 2
the vertex with id 4 connects to vertex with id 3
the vertex with id 5 connects to vertex with id 1
the vertex with id 6 connects to vertex with id 5

and this corresponds to a path

# HOW TO COMPILE:

g++ NetworkMotifCount.cpp -o NetworkMotifCount

change to long double if the number of trees is large but notice that it only works on 64 bit machine.

# HOW TO RUN:

./NetworkMotifCount 6trees 1 ExampleNetwork

the first parameter is the file containing all tree topoligies
the second parameter is the id of the tree we would like to count
the third parameter is the name of the network


