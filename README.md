# NETWORK MOTIF COUNTING

Network motif counting is a tool for counting the number of occurrences of a given undirected tree in a undirected network. The edges in the network could be have weights in [0,1]. It implements a randomized algorithm that are described in [our papers](#papers).

## COMPILATION

```
g++ NetworkMotifCount.cpp -o NetworkMotifCount
```

If the number of trees is large, you change the type to ```long double```. However, it only works on 64 bit machine.

## INPUT FILES
### Network File Format 

We include a sample network file ExampleNetwork. The first line contains the number of vertices and edges following by the edges together with their confidence score (should be between 0 and 1). All vertex ids start from 0.

### Tree Topology File Format 

Here we enclose all undirected unlabelled tree topologies with 6,7,8,9 and 10 vertices.

6trees containing all tree topologies of size 6 
7trees containing all tree topologies of size 7
...

The first number in each file is the number of vertices that all the trees in this file contain. Each line following describe a topology of a tree. The format of the tree topology is described through this following example:

0 1 2 3 1 5

There is an edge between vertex 1 and 2, vertex 2 and 3, vertex 3 and 4, vertex 1 and 5, vertex 5 and 6. This describes a simple path 6-5-1-2-3-4.

## USAGE

```
./NetworkMotifCount 6trees 1 ExampleNetwork
```

First parameter is the file containing all tree topoligies, followed by the id of the tree we would like to count and the file describing the network topology.

<a name="papers"></a>
## CITATIONS

Please cite the following papers:

1. Noga Alon, Phuong Dao, Iman Hajirasouliha, Fereydoun Hormozdiari, S. Cenk Sahinalp. "Biomolecular Network Motif Counting and Discovery by Color Coding". Bioinformatics (2008) 24 (13): i241-i249. [[LINK]](http://bioinformatics.oxfordjournals.org/content/24/13/i241.full)
2. Phuong Dao, Alexander Schonhuth, Fereydoun Hormozdiari, Iman Hajirasouliha, S. Cenk Sahinalp and Martin Ester. "Quantifying Systemic Evolutionary Changes by Color Coding Condence-Scored PPI Networks". Proceedings of the 9th Workshop on Algorithms on Bioinformatics (WABI 2009). [[LINK]](http://link.springer.com/chapter/10.1007%2F978-3-642-04241-6_4)

## CONTACTS

Please report any problems directly to the github issue tracker. Also, you can also send your feedbacks to phuongdao1@gmail.com.

## LICENSE

Network Motif Counting is distributed under GNU GPL license version 3.
