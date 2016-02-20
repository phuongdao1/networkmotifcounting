#include <cstdlib>
#include <sys/types.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string.h>

using namespace std;

int n;          // size of graph
int m;          // size of query graph
double eps;     // error probability
double delta;   // additive approximation constant
int numExp;     // the number of samples we need to take
double prob;
//int sameAsRoot=2; // the number of vertices that play the same role as root
double  overCounting=1.0;
int twoToPowerm;
int** subset;
//int mulFactor=5;
int* map;
bool* marked;



// struct for undirected graph
struct Graph
{
       int n;
       int** al;   // adjacent list
       char** am;   // adjacent matrix
       double** pm; // probability matrix
       bool** matchTo;
};

// struct for tree
// represented by a directed acylic graph with its root that has no node points to it
struct Tree
{
       int m;                 // the number of nodes
       int root;              // which node is the root of the tree
       int** al;              // adjacent list
       int** role;            // the number of siblings that a child may play the same role i.e. their rooted subtrees are isomorphic
       int* size;
       int* degree;
       //int** degreeSequence;
       int** ei;
       //int** ej;
};

struct Colorings
{
       int num;
       int** coloring;
       double * numColorings;
};


double * samplingTree;
int* samplingLoc;
double *** num;
Colorings c;
Graph g;
Tree t;
double RANDMAX=2147483648.0;
int range;
double first;
double second;


double ran()
{
       if (c.numColorings[c.num-1]<RANDMAX)
          return ((double)rand()*c.numColorings[c.num-1]/RANDMAX);
       else
       {
           range=(int)ceil(c.numColorings[c.num-1]/RANDMAX);
           second=((double)rand()/RANDMAX);
           first=(double)(rand()%range);
           if (first<range-1)
              return ((first+second)*RANDMAX);
           else
               return (first*RANDMAX+(c.numColorings[c.num-1]-first*RANDMAX)*second);
           //(double)(rand()%first)*(double)RAND_MAX
           //+second*c.numColorings[c.num-1]*(double)
       }
}

void swap(int *a, int *b)
{
  int t=*a; *a=*b; *b=t;
}

void swapd(double  *a, double  *b)
{
  double  t=*a; *a=*b; *b=t;
}

void quickSort(int a[],int lo, int hi)
{
    int i=lo, j=hi, h;
    int x=a[(lo+hi)/2];
    //double x1=a1[(lo+hi)/2];

    //  partition
    do
    {
        while (a[i]<x) i++;
        while (a[j]>x) j--;
        if (i<=j)
        {
            swap(&a[i],&a[j]);
            //swapd(&a1[i],&a1[j]);
            i++; j--;
        }
    } while (i<=j);

    //  recursion
    if (lo<j) quickSort(a, lo, j);
    if (i<hi) quickSort(a, i, hi);
}

void createSampleTree()
{
     //t.root=0;

     t.root=0;
     t.al[0][0]=2;
                  t.al[0][1]=1;
                  t.al[0][2]=4;
     t.al[1][0]=1;
                  t.al[1][1]=2;
     t.al[2][0]=1;
		    t.al[2][1]=3;
     t.al[4][0]=1;
                 t.al[4][1]=5;
     t.al[5][0]=1;
                  t.al[5][1]=6;
    /*
     t.al[0][0]=1;
                  t.al[0][1]=1;
     t.al[1][0]=1;
                  t.al[1][1]=2;
     t.al[2][0]=1;
                  t.al[2][1]=3;
     t.al[3][0]=1;
                 t.al[3][1]=4;
     t.al[4][0]=1;
                  t.al[4][1]=5;
     t.al[5][0]=1;
                  t.al[5][1]=6;
     t.al[6][0]=1;
                 t.al[6][1]=7;
     t.al[7][0]=1;
                  t.al[7][1]=8;
     t.al[8][0]=1;
                  t.al[8][1]=9;
     //t.al[9][0]=1;
     //             t.al[9][1]=10;
  */

     t.role[0][0]=1;
                    t.role[0][1]=1;
     t.role[1][0]=1;
                    t.role[1][1]=1;
     t.role[2][0]=1;
                    t.role[2][1]=1;
     t.role[3][0]=1;
                    t.role[3][1]=1;
     //t.role[4][0]=1;
     //               t.role[4][1]=1;
     //t.role[5][0]=1;
     //               t.role[5][1]=1;
     //t.role[6][0]=0;
    // t.role[6][0]=1;
    //                t.role[6][1]=1;
    // t.role[7][0]=1;
    //                t.role[7][1]=1;
   //  t.role[8][0]=0;
                    //t.role[8][1]=1;
     /*
     t.al[2][0]=3;
                  t.al[2][1]=1;
                  t.al[2][2]=0;
                  t.al[2][3]=3;
     t.al[0][0]=1;
                  t.al[0][1]=4;

     t.role[2][0]=3;
                  t.role[2][1]=1;
                  t.role[2][2]=1;
                  t.role[2][3]=2;
     t.role[0][0]=1;
                  t.role[0][1]=1;
     */

     /*
     t.al[2][0]=3;
                  t.al[2][1]=1;
                  t.al[2][2]=0;
                  t.al[2][3]=3;
     t.al[0][0]=2;
                  t.al[0][1]=4;
                  t.al[0][2]=5;


     t.role[2][0]=3;
                  t.role[2][1]=1;
                  t.role[2][2]=1;
                  t.role[2][3]=2;
     t.role[0][0]=2;
                  t.role[0][1]=1;
                  t.role[0][2]=2;
     */

     /*
     t.al[2][0]=3;
                  t.al[2][1]=1;
                  t.al[2][2]=0;
                  t.al[2][3]=3;
     t.al[0][0]=2;
                  t.al[0][1]=4;
                  t.al[0][2]=5;
     t.al[8][0]=3;
                  t.al[8][1]=7;
                  t.al[8][2]=6;
                  t.al[8][3]=9;
     t.al[6][0]=2;
                  t.al[6][1]=10;
                  t.al[6][2]=11;
     t.al[12][0]=2;
                   t.al[12][1]=2;
                   t.al[12][2]=8;


     t.role[2][0]=3;
                  t.role[2][1]=1;
                  t.role[2][2]=1;
                  t.role[2][3]=2;
     t.role[0][0]=2;
                  t.role[0][1]=1;
                  t.role[0][2]=2;
     t.role[8][0]=3;
                  t.role[8][1]=1;
                  t.role[8][2]=1;
                  t.role[8][3]=2;
     t.role[6][0]=2;
                  t.role[6][1]=1;
                  t.role[6][2]=2;
     t.role[12][0]=2;
                  t.role[12][1]=1;
                  t.role[12][2]=2;
     */

     /*
     t.al[2][0]=4;
                  t.al[2][1]=1;
                  t.al[2][2]=0;
                  t.al[2][3]=3;
                  t.al[2][4]=4;


     t.role[2][0]=4;
                  t.role[2][1]=1;
                  t.role[2][2]=2;
                  t.role[2][3]=3;
                  t.role[2][4]=4;
     */
     /*
     t.al[2][0]=3;
                  t.al[2][1]=1;
                  t.al[2][2]=0;
                  t.al[2][3]=3;



     t.role[2][0]=3;
                  t.role[2][1]=1;
                  t.role[2][2]=2;
                  t.role[2][3]=3;
     */

}

void createGraph(Graph* g, int n)
{
     int i,j;
     g->n=n;
     g->al=(int**)malloc(sizeof(int*)*n);
      g->matchTo=(bool**)malloc(sizeof(bool*)*n);

     //g->am=(char**)malloc(sizeof(char*)*n);
     g->pm=(double**)malloc(sizeof(double*)*n);
     for (i=0;i<n;i++)
     {
         g->al[i]=(int*)malloc(sizeof(int)*(n+1));
         g->matchTo[i]=(bool*)malloc(sizeof(bool)*(m));
         //g->am[i]=(char*)malloc(sizeof(char)*n);
         g->pm[i]=(double*)malloc(sizeof(double)*n);

         g->al[i][0]=0;
     }
     for (i=0;i<n;i++)
     for (j=0;j<m;j++)
         g->matchTo[i][j]=false;

     for (i=0;i<n;i++)
     for (j=0;j<n;j++)
         g->pm[i][j]=0.0;
}

void deleteGraph(Graph* g, int n)
{
     int i;
     for (i=0;i<n;i++)
     {
         free(g->al[i]);
         free(g->matchTo[i]);
         //free(g->am[i]);
         free(g->pm[i]);
     }
     free(g->al);
     free(g->matchTo);
     free(g->pm);
     //free(g->am);
     //free(g->degreeSequence);
}

void addEdge(Graph* g,int u,int v,double prob)
{
     //cout<<g->al[u][0]<<" "<<g->al[v][0]<<endl;
     g->al[u][0]++;
     g->al[v][0]++;
     g->al[u][ g->al[u][0] ]=v;
     g->al[v][ g->al[v][0] ]=u;
     //g->al[u][0]++;
     //g->al[v][0]++;
     //g->am[u][v]=1;
     //g->am[v][u]=1;
     g->pm[u][v]=prob;
     g->pm[v][u]=prob;

}

void addEdge1(Graph* g,int u,int v)
{
     //cout<<g->al[u][0]<<" "<<g->al[v][0]<<endl;
     g->al[u][0]++;
     g->al[v][0]++;
     g->al[u][ g->al[u][0] ]=v;
     g->al[v][ g->al[v][0] ]=u;
     //g->al[u][0]++;
     //g->al[v][0]++;
     //g->am[u][v]=1;
     //g->am[v][u]=1;
     //g->pm[u][v]=prob;
     //g->pm[v][u]=prob;

}

void addTreeEdge(Tree* g,int u,int v)
{
     g->al[u][ g->al[u][0] ]=v;
     g->al[u][0]++;
}

void createTree(Tree* g, int m)
{
     int i;
     g->m=m;
     g->al=(int**)malloc(sizeof(int*)*m);
     //g->degreeSequence=(int**)malloc(sizeof(int*)*m);
     g->role=(int**)malloc(sizeof(int*)*m);
     g->degree=(int*)malloc(sizeof(int)*m);
     g->size=(int*)malloc(sizeof(int)*m);
     g->ei=(int**)malloc(sizeof(int*)*m);
     //g->ej=(int**)malloc(sizeof(int*)*m);
     for (i=0;i<m;i++)
     {
         g->al[i]=(int*)malloc(sizeof(int)*m);
         //g->degree[i]=(int*)malloc(sizeof(int)*m);
         g->role[i]=(int*)malloc(sizeof(int)*m);
         g->ei[i]=(int*)malloc(sizeof(int)*(n+1));
         //g->ej[i]=(int*)malloc(sizeof(int)*n*mulFactor);
         g->ei[i][0]=0;
         //g->ej[i][0]=0;
         //g->degreeSequence[i]=(int*)malloc(sizeof(int)*(m+1));
         //g->degreeSequence[i][0]=0;
         //g->size[i]=(int*)malloc(sizeof(int)*m);
         g->role[i][0]=0;
         g->al[i][0]=0;
     }
}

void deleteTree(Tree* g, int m)
{
     int i;
     for (i=0;i<m;i++)
     {
         free(g->al[i]);
        // cout<<i<<" normal "<<m<<endl;
         free(g->role[i]);
        // cout<<i<<" normal "<<m<<endl;
         //free(g->degreeSequence[i]);
         //cout<<i<<" normal "<<m<<endl;
         free(g->ei[i]);
        // cout<<i<<" normal "<<m<<endl;
         //free(g->ej[i]);
         //free(g->degree[i]);
         //free(g->size[i]);
     }
    // cout<<"error here"<<endl;
     free(g->al);
     free(g->role);
     free(g->ei);
   //  cout<<"error here"<<endl;
     //free(g->ej);
     free(g->degree);
     //free(g->degreeSequence);
     free(g->size);
}

int computeSize(int root)
{
        t.size[root]=0;
        for (int i=1;i<=t.al[root][0];i++)
            t.size[root]+=computeSize(t.al[root][i]);
        t.size[root]++;
        return t.size[root];
}

void computeDegree(int root)
{
    if (root==t.root)
    {
       t.degree[root]=t.al[root][0];
       for (int i=1;i<=t.al[root][0];i++)
           computeDegree(t.al[root][i]);
    }
    else
    {
       t.degree[root]=t.al[root][0]+1;
       for (int i=1;i<=t.al[root][0];i++)
           computeDegree(t.al[root][i]);
    }
}

/*
void computeEi(int root)
{
     int i,k,l,i1,j1,k1;
     for (i=0;i<n;i++)
     if (t.degree[root]<=g.al[i][0])
     {
             t.ei[root][0]++;
             t.ei[root][t.ei[root][0]]=i;
     }
     for (i=1;i<=t.al[root][0];i++)
         computeEi(t.al[root][i]);
}
*/


bool match(int* a,int* b)
{
     if (a[0]==0)
        return true;
     if (b[0]==0)
        return false;
     int i,j;

     j=0;
     for (i=1;i<=a[0];i++)
     {
         do
         {
               j=j+1;
         } while ((j+1<=b[0])&&(a[i]>b[j]));
         if (j>b[0])
            return false;
         if (a[i]>b[j])
            return false;
     }
     return true;
}

void printSeq(int* arr)
{
     int i;
     cout<<"*";
     for (i=1;i<=arr[0];i++)
         cout<<arr[i]<<" ";
     cout<<"*"<<endl;
}

void computeEi(int root,int parent)
{
     int *a;
     int *b;
     a=(int*)malloc(sizeof(int)*m);
     b=(int*)malloc(sizeof(int)*n);
     int i,j,k,l,i1,j1,k1;
     t.ei[root][0]=0;
     for (i=0;i<n;i++)
     if (t.degree[root]<=g.al[i][0])
     {
             a[0]=0;
             b[0]=0;
             if (root!=t.root)
             {
                a[0]=1;
                a[1]=t.degree[parent];
             }
             for (j=1;j<=t.al[root][0];j++)
             {
                 a[0]=a[0]+1;
                 a[a[0]]=t.degree[t.al[root][j]];
             }
             for (j=1;j<=g.al[i][0];j++)
             {
                 b[0]++;
                 b[b[0]]=g.al[g.al[i][j]][0];
             }
             //printSeq(a);
             //printSeq(b);
             quickSort(a,1,a[0]);
             quickSort(b,1,b[0]);
             //printSeq(a);
             //printSeq(b);
             if (match(a,b))
             {
               t.ei[root][0]++;
               t.ei[root][t.ei[root][0]]=i;
               g.matchTo[i][root]=true;
             }
     }
     free(a);
     free(b);
     /*
     cout<<root<<" match "<<endl;
     for (i=1;i<=t.ei[root][0];i++)
         cout<<t.ei[root][i]<<"  ";
     cout<<endl;*/
     for (i=1;i<=t.al[root][0];i++)
         computeEi(t.al[root][i],root);
}

/*
int computeDegreeSequence(int root)
{
    if (root==t.root)
    {
       t.degree[root]=t.al[root][0];
       for (int i=1;i<=t.al[root][0];i++)
           computeDegree(t.al[root][i]);
    }
    else
    {
       t.degree[root]=t.al[root][0]+1;
       for (int i=1;i<=t.al[root][0];i++)
           computeDegree(t.al[root][i]);
    }
}
*/





void generateColorings(Colorings* c,int n,int m)
{
     int i,j,temp;
     int* arr;
     //temp=(int)ceil(pow(2.0,(m*1.0))/2.7);
     //temp=(int)ceil(pow(2.0,(m*1.0)/2.0));
     //c->num=(int)ceil(temp*(log(n*1.0)/log(2.0)));
     c->num=numExp;

     //cout<<"num of colorings: "<<c->num<<endl;
     //c->num=10;
     c->coloring=(int**)malloc(sizeof(int*)*c->num);
     c->numColorings=(double *)malloc(sizeof(double )*c->num);
     arr=(int*)malloc(sizeof(int)*m);
     int i1,j1,same;
     //srand((unsigned)time(0));
     //srand(1990);
     //srand(1985); //star 3 leaves
     for (i=0;i<c->num;i++)
     {
         c->coloring[i]=(int*)malloc(sizeof(int)*n);
         for (j=0;j<n;j++)
         {
             c->coloring[i][j]=rand()%m;
            // cout<<c->coloring[i][j];
         }
	 /*
         for (i1=0;i1<m;i1++)
         {
             do
             {
               arr[i1]=rand()%n;
               same=0;
               for (j1=0;j1<i1;j1++)
               if (arr[j1]==arr[i1])
               {
                 same=1; break ;
               }
             } while (same==1);
             c->coloring[i][arr[i1]]=i1;
         }
	 */
         /*
         cout<<"c "<<i<<" ";
         for (j=0;j<n;j++)
         {
            // c->coloring[i][j]=rand()%m;
             cout<<c->coloring[i][j];
         }

        cout<<endl;
        */
     }
     free(arr);

}


void reGenerateColorings(Colorings* c,int n,int m)
{
     int i,j,temp;
     int* arr;
     //temp=(int)ceil(pow(2.0,(m*1.0))/2.7);
     //temp=(int)ceil(pow(2.0,(m*1.0)/2.0));
     //c->num=(int)ceil(temp*(log(n*1.0)/log(2.0)));
     //cout<<"num of colorings: "<<c->num<<endl;
     //c->num=10;
     //c->coloring=(int**)malloc(sizeof(int*)*c->num);
     //c->numColorings=(double *)malloc(sizeof(double )*c->num);
     arr=(int*)malloc(sizeof(int)*m);
     int i1,j1,same;
     //srand((unsigned)time(0));
     //srand(1990);
     //srand(1985); //star 3 leaves
     for (i=0;i<c->num;i++)
     {
         //c->coloring[i]=(int*)malloc(sizeof(int)*n);
         for (j=0;j<n;j++)
         {
             c->coloring[i][j]=rand()%m;
            // cout<<c->coloring[i][j];
         }
	 /*
         for (i1=0;i1<m;i1++)
         {
             do
             {
               arr[i1]=rand()%n;
               same=0;
               for (j1=0;j1<i1;j1++)
               if (arr[j1]==arr[i1])
               {
                 same=1; break ;
               }
             } while (same==1);
             c->coloring[i][arr[i1]]=i1;
         }
	 */
         /*
         cout<<"c "<<i<<" ";
	 
         for (j=0;j<n;j++)
         {
            // c->coloring[i][j]=rand()%m;
             cout<<c->coloring[i][j];
         }

        cout<<endl;
	 */
        
     }




     free(arr);

}

void deleteColorings(Colorings* c,int n,int m)
{
     int i;
     for (i=0;i<c->num;i++)
         free(c->coloring[i]);
     free(c->coloring);
     free(c->numColorings);
}

#define ranf() ((double)(rand()%32768+1)/32768.0)

double  countTrees(/*Colorings* c,*/int id,int root,int size)
{
     //memset(num,0,sizeof(int)*n*m*twoToPowerm);
     //cout<<"    "<<root<<" "<<size<<endl;
     double  total=0.0;
     if (t.al[root][0]==0)
     {
        //cout<<"cccc"<<root<<" "<<size<<endl;
        /*
        for (int i=0;i<n;i++)
        if (t.degree[root]<=g.al[i][0])
        {
            //cout<<"iii"<<(1<<c->coloring[id][i])<<endl;
            num[i][root][1<<c.coloring[id][i]]++;
            //cout<<"iii"<<i<<endl;
            total++;
        }
        */
//cout<<"don be here"<<endl;
        for (int i=1;i<=t.ei[root][0];i++)
        {
            num[t.ei[root][i]][root][1<<c.coloring[id][t.ei[root][i]]]++;
            total++;
        }
        //return total;
     }
     else
     {
         //cout<<"dddd"<<root<<" "<<size<<endl;
  //       cout<<"overflow here111"<<endl;
	 int mem=t.al[root][t.al[root][0]];
         int subtreeSize=t.size[mem];
         int sameRole=t.role[root][t.role[root][0]];
         int i,j,k,test,anotherSet;
         int u,v;
         //cout<<"xxx"<<mem<<" "<<subtreeSize<<" "<<sameRole<<endl;
       //  cout<<"overflow here"<<endl;
	 t.al[root][t.al[root][0]]=0;
         t.al[root][0]--;
         t.role[root][0]--;


         countTrees(/*c,*/id,root,size-subtreeSize);
         countTrees(/*c,*/id,mem,subtreeSize);

         t.role[root][0]++;
         t.al[root][0]++;
         t.al[root][t.al[root][0]]=mem;
         int uu;
         int vv;

         //if ((subset[size][0]*subset[subtreeSize][0]<subset[subtreeSize][0]*subset[size-subtreeSize][0])&&(subset[size][0]*subset[subtreeSize][0]<subset[size][0]*subset[size-subtreeSize][0]))
          //cout<<"overflow here222"<<" "<<size<<" "<<subtreeSize<<endl;
         // int stupid=-subtreeSize+size;
          //cout<<subset[stupid][size-subtreeSize]<<" "<<size<<" "<<subtreeSize<<endl;

         //if ((subset[size][0]<subset[size-subtreeSize][0])&&(subset[subtreeSize][0]<subset[size-subtreeSize][0]))
         if (subset[size][0]<subset[size-subtreeSize][0])/*&&(subset[subtreeSize][0]<subset[size-subtreeSize][0]))*/
         {
	///	 cout<<"overflow herexxx"<<endl;

            for (uu=1;uu<=t.ei[root][0];uu++)
            {
                u=t.ei[root][uu];
                //if (t.degree[root]<=g.al[u][0])
                for (i=1;i<=subset[size][0];i++)
                {
                    for (k=1;k<=g.al[u][0];k++)
                    if (c.coloring[id][u]!=c.coloring[id][g.al[u][k]])
                    {
                       int v=g.al[u][k];
                       //if (t.degree[mem]<=g.al[v][0])                 //yes yes
                       if (g.matchTo[v][mem])
                       for (j=1;j<=subset[subtreeSize][0];j++)
                       {
                           anotherSet=subset[size][i]^subset[subtreeSize][j];
                           test=subset[size][i]|subset[subtreeSize][j];
                           if (test==subset[size][i])
                       //{
                              num[u][root][subset[size][i]]+=(num[v][mem][subset[subtreeSize][j]]*num[u][root][anotherSet])*g.pm[u][v];
                          //if ((num[v][mem][subset[subtreeSize][j]]>0)&&(num[u][root][anotherSet]>0))
                            // cout<<"mmm"<<u<<" "<<root<<" "<<v<<" "<<mem<<" "<<" "<<num[u][root][anotherSet]<<" "<<num[v][mem][subset[subtreeSize][j]]<<endl;
                       //}
                       }
                    }
                //num[u][root][subset[size][i]]/=sameRole;
                    total+=num[u][root][subset[size][i]];
                }
            }
         }
         else /*if (subset[subtreeSize][0]<subset[size][0])*/
         {
       //  cout<<"go here sometimes"<<endl;
          for (uu=1;uu<=t.ei[root][0];uu++)
          {
            u=t.ei[root][uu];
            //if (t.degree[root]<=g.al[u][0])
            for (i=1;i<=subset[size-subtreeSize][0];i++)
            if (num[u][root][subset[size-subtreeSize][i]]>0)
            {
                for (k=1;k<=g.al[u][0];k++)
                if (c.coloring[id][u]!=c.coloring[id][g.al[u][k]])
                {
                    int v=g.al[u][k];
                    //if (t.degree[mem]<=g.al[v][0])                //yes yes
                    if (g.matchTo[v][mem])
                    for (j=1;j<=subset[subtreeSize][0];j++)
                    if (num[v][mem][subset[subtreeSize][j]]>0)
                    {
                       anotherSet=subset[size-subtreeSize][i]|subset[subtreeSize][j];
                       test=subset[size-subtreeSize][i]&subset[subtreeSize][j];
                       if (test==0)
                       //{
                       {
                          num[u][root][anotherSet]+=(num[v][mem][subset[subtreeSize][j]]*num[u][root][subset[size-subtreeSize][i]])*g.pm[u][v];

                          total+=(num[v][mem][subset[subtreeSize][j]]*num[u][root][subset[size-subtreeSize][i]])*g.pm[u][v];
                       }
                          //if ((num[v][mem][subset[subtreeSize][j]]>0)&&(num[u][root][anotherSet]>0))
                            // cout<<"mmm"<<u<<" "<<root<<" "<<v<<" "<<mem<<" "<<" "<<num[u][root][anotherSet]<<" "<<num[v][mem][subset[subtreeSize][j]]<<endl;
                       //}
                    }
                }
                //num[u][root][subset[size][i]]/=sameRole;
                //total+=num[u][root][subset[anotherSet][i]];
              }
            }
           }
         /*else
         {
         // cout<<"it is clear here"<<endl;
          for (uu=1;uu<=t.ei[root][0];uu++)
          {
            u=t.ei[root][uu];
            //if (t.degree[root]<=g.al[u][0])
            for (i=1;i<=subset[size-subtreeSize][0];i++)
            if (num[u][root][subset[size-subtreeSize][i]]>0)
            {
                for (k=1;k<=g.al[u][0];k++)
                if (c.coloring[id][u]!=c.coloring[id][g.al[u][k]])
                {
                    int v=g.al[u][k];
                    //if (t.degree[mem]<=g.al[v][0])                //yes yes
                    if (g.matchTo[v][mem])
                    for (j=1;j<=subset[size][0];j++)
                    {
                       anotherSet=subset[size][j]^subset[size-subtreeSize][i];
                       test=subset[size-subtreeSize][i]|subset[size][j];
                       if (test==subset[size][j])
                       //{
                       {
                          num[u][root][subset[size][j]]+=(num[v][mem][anotherSet]*num[u][root][subset[size-subtreeSize][i]]);
                          total+=(num[v][mem][anotherSet]*num[u][root][subset[size-subtreeSize][i]]);
                       }
                          //if ((num[v][mem][subset[subtreeSize][j]]>0)&&(num[u][root][anotherSet]>0))
                            // cout<<"mmm"<<u<<" "<<root<<" "<<v<<" "<<mem<<" "<<" "<<num[u][root][anotherSet]<<" "<<num[v][mem][subset[subtreeSize][j]]<<endl;
                       //}
                    }
                }
                //num[u][root][subset[size][i]]/=sameRole;
                //total+=num[u][root][subset[anotherSet][i]];
            }
         }



         }
         */
     }
     //if (total>2000)
     //   cout<<"YES YES YES"<<endl;
     //cout<<" "<<root<<" "<<size<<" "<<total<<endl;
     return total;
}

void sort(int a[],double  a1[], int lo, int hi)
{
    int i=lo, j=hi, h;
    int x=a[(lo+hi)/2];
    double  x1=a1[(lo+hi)/2];

    //  partition
    do
    {
        while ((a[i]<x)||((a[i]==x)&&(a1[i]<x1))) i++;
        while ((a[j]>x)||((a[j]==x)&&(a1[j]>x1))) j--;
        if (i<=j)
        {
            swap(&a[i],&a[j]);
            swapd(&a1[i],&a1[j]);
            i++; j--;
        }
    } while (i<=j);

    //  recursion
    if (lo<j) sort(a, a1, lo, j);
    if (i<hi) sort(a, a1, i, hi);
}

void sortMidArr(double  a[], int lo, int hi)
{
    int i=lo, j=hi, h;
    //int x=a[(lo+hi)/2];
    double  x=a[(lo+hi)/2];

    //  partition
    do
    {
        while (a[i]<x) i++;
        while (a[j]>x) j--;
        if (i<=j)
        {
            //swap(&a[i],&a[j]);
            swapd(&a[i],&a[j]);
            i++; j--;
        }
    } while (i<=j);

    //  recursion
    if (lo<j) sortMidArr(a, lo, j);
    if (i<hi) sortMidArr(a, i, hi);
}

double countAllTrees(/*Colorings* c*/)
{
    int i,j,k,l;
    //int midExp=(int)ceil(log(1.0/delta));
    int midExp=(int)ceil(log(1/delta)); // the value of log (1/delta)
    double midArr[100];

    //int midExp=log(1.0/delta);

    for (l=0;l<midExp;l++)
    {
        //c.numColorings[-1]=0;
        //cout<<"mid "<<l<<endl;
      reGenerateColorings(&c,n,m);
      for (k=0;k<c.num;k++)
        {
            //cout<<"num "<<k<<endl;
            for (i=0;i<n;i++)
            for (j=0;j<m;j++)
                memset(num[i][j],0,sizeof(double )*twoToPowerm);

            //cout<<"counting "<<k<<"..."<<endl;
            //cout<<"before "<<t.root<<" "<<t.size[t.root]<<" "<<m<<endl;
            c.numColorings[k]=countTrees(/*c,*/k,t.root,t.size[t.root]);
	    int jj;
           for (jj=0;jj<n;jj++)
         {
	   //printf("%lf %lf",c.num,numExp);
            // c->coloring[i][j]=rand()%m;
	   //cout<<c.coloring[k][jj];
         }

	   //cout<<endl;
	   //printf(" %lf\n",c.numColorings[k]);
            //cout<<"num "<<k<<" "<<setprecision(30)<<c.numColorings[k]<<endl;
            if (k>0)
            c.numColorings[k]+=c.numColorings[k-1];
            //cout<<"num trees "<<countTrees(&c,k,t.root,t.size[t.root])<<endl;
        }

      midArr[l]=(c.numColorings[c.num-1]/(c.num*1.0))*(1/prob);
      //printf("%lf\n",midArr[l]);
        //cout<<"mid "<<l<<endl;
        //cout<<setprecision(30)<<"Total trees: "<<midArr[l]<<endl;
    }

    //cout<<setprecision(30)<<"Total trees: "<<midArr[midExp/2]<<endl;
    //cout<<"Total trees: "<<midArr[c.num/2];
    //printf("Total trees: %lf\n",midArr[c.num/2]);
    sortMidArr(midArr,0,midExp-1);
    return midArr[midExp/2];
     //for (k=0;k<c.num;k++)
     //    cout<<"num "<<k<<" "<<setprecision(20)<<c.numColorings[k]<<endl;
}


void readFile()
{
    ifstream inFile;



    inFile.open("Ecoli");
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    int n,numEdges;
    inFile>>n>>numEdges;

//    while (inFile>>x>>y)
        //addEdge(&g,1,0);

    inFile.close();
}

Graph gt;

double overCount(int* c,int root,int size,Graph* g)
{
     double total=0;
     if (t.al[root][0]==0)
     {
        //cout<<"going here"<<endl;
        for (int i=0;i<m;i++)
        if (t.degree[root]<=g->al[i][0])
        {
            num[i][root][1<<c[i]]++;
            total++;
        }
        //cout<<"error going here"<<endl;
     }
     else
     {
         //cout<<root<<" "<<size<<endl;
         int mem=t.al[root][t.al[root][0]];
         int subtreeSize=t.size[mem];
        // int sameRole=t.role[root][t.role[root][0]];
         int i,j,k,test,anotherSet;
         int u,v;

         t.al[root][t.al[root][0]]=0;
         t.al[root][0]--;
         //t.role[root][0]--;

         overCount(c,root,size-subtreeSize,g);
         overCount(c,mem,subtreeSize,g);

         //t.role[root][0]++;
         t.al[root][0]++;
         t.al[root][t.al[root][0]]=mem;
         int uu;
         int vv;


            for (uu=0;uu<m;uu++)
            {
                u=uu;
                if (t.degree[root]<=g->al[u][0])
                for (i=1;i<=subset[size][0];i++)
                {
                    for (k=1;k<=g->al[u][0];k++)
                    if (c[u]!=c[g->al[u][k]])
                    {
                       int v=g->al[u][k];
                       if (t.degree[mem]<=g->al[v][0])
                       for (j=1;j<=subset[subtreeSize][0];j++)
                       {
                           anotherSet=subset[size][i]^subset[subtreeSize][j];
                           test=subset[size][i]|subset[subtreeSize][j];
                           if (test==subset[size][i])

                              num[u][root][subset[size][i]]+=(num[v][mem][subset[subtreeSize][j]]*num[u][root][anotherSet]);

                       }
                    }

                    total+=num[u][root][subset[size][i]];
                }

         }
     }
     return total;
}

void copyTree(int root)
{
    for (int i=1;i<=t.al[root][0];i++)
    {
           //cout<<"Add edge "<<root<<" "<<t.al[root][i]<<endl;
           addEdge1(&gt,root,t.al[root][i]);
           copyTree(t.al[root][i]);
    }
}


void computeOvercountingFactor()
{
       int Coloring[30];
       int i,j;
       for (i=0;i<m;i++)
           Coloring[i]=i;
       createGraph(&gt,m);
       copyTree(t.root);
       //cout<<"Error up there."<<endl;
        for (i=0;i<m;i++)
        for (j=0;j<m;j++)
            memset(num[i][j],0,sizeof(double)*twoToPowerm);

       overCounting=overCount(Coloring,t.root,m,&gt)*1.0;
       //cout<<"Overcounting factor: "<<overCounting<<endl;
       deleteGraph(&gt,m);
}

#define delims " \t"

int main(int argc, char *argv[])
{
    ifstream inFile;


    /*    for (ii=0;ii<m;ii++)
    {
        for (jj=1;jj<=t.al[ii][0];jj++)
            cout<<" "<<t.al[ii][jj];
        cout<<endl;
    }*/
    inFile.open(argv[3]);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    int numEdges;
    //m=7;

    inFile>>n>>numEdges;
    //n=300;
    //cout<<n<<" "<<numEdges<<endl;
    createGraph(&g,n);

    //while (inFile >> x>>" ">>y)
    for (int ttt=0;ttt<numEdges;ttt++)
    {
        int x,y;
        double prob;
        inFile >>x>>y>>prob;
        //cout<<x<<" "<<y<<endl;
        addEdge(&g,x,y,prob);
    }
    inFile.close();

    //int ii,jj;


    FILE *f;
    FILE *f1;
    char line[1024];
    char resultFile[256];
    int treeId;
    int countId=0;
    int ii,jj;
    char *token;
    int numComp;
    sscanf (argv[2], "%d", &treeId);



    f=fopen(argv[1],"rt");

    fgets(line, 1024, f);
    sscanf (line, "%d", &m);
    sprintf(resultFile,"%s_%d_%d",argv[3],m,treeId);
    twoToPowerm=(int)pow(2.0,m*1.0);
    delta=0.01;
    eps=0.1;
    createTree(&t,m);

    //cout<<"xxx "<<m<<" "<<treeId<<endl;
    do
    {
      fgets(line, 1024, f);
      //cout<<line<<endl;
      //cout<<line<<" "<<countId+1<<endl;
      countId++;
    } while (countId<treeId);
      token=strtok( line, delims );
      token=strtok( NULL, delims );
      ii=0;
      while (token!=NULL)
      {
            ii++;
            sscanf (token, "%d", &jj);
            jj--;
            t.al[jj][0]++;
            t.al[jj][t.al[jj][0]]=ii;
            //cout<<"eee "<<jj<<" "<<ii<<endl;
            if (ii==m-1)
		break;
            token=strtok( NULL, delims );
      }
      t.root=0;
    //createSampleTree();
    fclose(f);

    /*
    cout<<"list edges"<<endl;
        for (ii=0;ii<m;ii++)
    {
        for (jj=1;jj<=t.al[ii][0];jj++)
            cout<<" "<<t.al[ii][jj];
        cout<<endl;
    }
   */
    //cout<<"error down here"<<endl;
   /*
    for (int x=0;x<n;x++)
    for (int y=x+1;y<n;y++)
        addEdge(&g,x,y);
    */
    //n=13;
    //m=6;
    //twoToPowerm=(int)pow(2.0,m*1.0);
    //delta=0.1;
    //eps=0.1;

    srand((unsigned)time(0));

    //createGraph(&g,n);
    //createTree(&t,m);
    int i,j,k;
    prob=1.0;
    for (i=1;i<=m;i++)
      prob=prob*(i*1.0);
    for (i=1;i<=m;i++)
      prob=prob/(m*1.0);
    numExp=(int)ceil(4*(1.0/prob)/(eps*eps));
    generateColorings(&c,n,m);
    //int i,j,k;

    //cout<<"error here ???"<<endl;

    num=(double ***)malloc(sizeof(double **)*n);
    map=(int*)malloc(sizeof(int)*m);
    marked=(bool*)malloc(sizeof(bool)*m);
    for (i=0;i<n;i++)
        num[i]=(double **)malloc(sizeof(double *)*m);
    for (i=0;i<n;i++)
        for (j=0;j<m;j++)
            num[i][j]=(double *)malloc(sizeof(double )*twoToPowerm);

    subset=(int**)malloc(sizeof(int*)*(m+1));
    for (i=0;i<=m;i++)
    {
        subset[i]=(int*)malloc(sizeof(int)*twoToPowerm);
        subset[i][0]=0;
    }
    int test;
    int count1;

    for (int i=1;i<twoToPowerm;i++)
    {
         count1=0;
         for (j=0;j<m;j++)
         {
             test = i & (1 << j);
             if (test>0)
                count1++;
         }
         subset[count1][0]++;
         subset[count1][subset[count1][0]]=i;
    }

    //cout<<"error here ???"<<endl;
    //cout<<"get here"<<endl;
    /*
    for (i=7;i>=5;i--)
        for (j=i-1;j>=4;j--)
        {
            //cout<<i<<" "<<j<<endl;
            addEdge(&g,i,j);
        }
    for (i=3;i>=1;i--)
        for (j=i-1;j>=0;j--)
        {
            //cout<<i<<" "<<j<<endl;
            addEdge(&g,i,j);
        }
    addEdge(&g,4,0);
    */

    /*
    for (i=12;i>=10;i--)
        for (j=i-1;j>=9;j--)
            addEdge(&g,i,j);

    for (i=8;i>=6;i--)
        for (j=i-1;j>=5;j--)
            addEdge(&g,i,j);

    for (i=4;i>=2;i--)
    for (j=i-1;j>=1;j--)
            addEdge(&g,i,j);


    addEdge(&g,1,0);
    addEdge(&g,5,0);
    addEdge(&g,9,0);
    */


    //cout<<"error here ???"<<endl;
    /*
    bool** edge;
    edge=(bool**)malloc(sizeof(bool*)*n);
    for (i=0;i<n;i++)
        edge[i]=(bool*)malloc(sizeof(bool)*n);

    //bool edge[n][n];
    //cout<<"error here ???"<<endl;


    for (i=0;i<n;i++)
        for (j=0;j<n;j++)
            edge[i][j]=false;
    for (k=1;k<=4*n;k++)
    {
        do
        {
        do
        {
           i=rand()%n;
           j=rand()%n;
        } while (i==j);
        } while (edge[i][j]);
        edge[i][j]=true;
        addEdge(&g,i,j);
    }
    */

    //cout<<"error here ???"<<endl;


    //numExp=(int)ceil((c.num*4.0)*log(2.0/delta)/(eps*eps));

    samplingTree=(double *)malloc(sizeof(double )*numExp);
    samplingLoc=(int*)malloc(sizeof(int)*numExp);

    //cout<<"error here ???"<<endl;
    //createSampleTree();
    computeSize(t.root);
    computeDegree(t.root);
    //cout<<"error here ???"<<endl;
    //computeEi(t.root);

    computeEi(t.root,-1);
    //cout<<"error some where here"<<endl;
    // cout<<"up here"<<endl;
    /*
    for (i=0;i<n;i++)
        for (j=0;j<m;j++)
            memset(num[i][j],0,sizeof(int)*twoToPowerm);
    //cout<<"taman"<<endl;
    cout<<"num trees "<<countTrees(&c,1,t.root,t.size[t.root])<<endl;

    for (i=0;i<n;i++)
        for (j=0;j<m;j++)
            memset(num[i][j],0,sizeof(int)*twoToPowerm);
    //cout<<"taman"<<endl;
    cout<<"num trees "<<countTrees(&c,1,t.root,t.size[t.root])<<endl;
    */
    //countAllTrees(&c);
    /*
    for (i=0;i<n;i++)
    for (j=0;j<m;j++)
        memset(num[i][j],0,sizeof(int)*twoToPowerm);

    int numTrees=countTrees(&c,1,t.root,t.size[t.root]);
    //int map[m];
    cout<<"Num of trees: "<<numTrees<<endl;
    for (i=0;i<numTrees;i++)
    {
        for (j=0;j<m;j++)
            map[j]=-1;
        //sampleTrees1(&c,1,t.root,m,i+1,map,twoToPowerm-1);
        int total=0,pre=0;
        for (int u=0;u<n;u++)
        {
            pre=total;
            total+=num[u][t.root][twoToPowerm-1];
            if (total>=i+1)
            {
               sampleTrees2(&c,1,t.root,u,m,i+1-pre,map,twoToPowerm-1);
               break;
            }
        }


        for (j=0;j<m;j++)
            cout<<" "<<map[j];
        cout<<endl;
    }
    cout<<"test "<<num[4][0][24]<<endl;
     */

    /*
    for (i=0;i<m;i++)
    {
        for (j=1;j<=t.ei[i][0];j++)
            cout<<"e "<<i<<" "<<t.ei[i][j]<<" "<<t.ej[i][j]<<endl;
    }*/
    /*
    for (i=0;i<m;i++)
        cout<<"d "<<i<<" "<<t.degree[i]<<endl;
    for (i=0;i<m;i++)
        cout<<"s "<<i<<" "<<t.size[i]<<endl;*/

   computeOvercountingFactor();
   //cout<<"counting here..."<<endl;
   int start=(unsigned)time(0);
   srand((unsigned)time(0));
   double appResult=countAllTrees()/(overCounting*1.0);
   //cout<<"Finishing counting trees"<<endl;
   //cout<<"# Occurrences: "<<setprecision(30)<<sampling(&g)<<endl;
   //double appResult=sampling(&g);
   //cout<<"# Occurrences: "<<setprecision(30)<<appResult<<endl;
   int ending=(unsigned)time(0);
   //cout<<-start+ending<<" seconds."<<endl;

   //f1=fopen(resultFile,"a");
   //fprintf(f1,"%d %.0lf %.0lf %d\n",treeId,appResult,overCounting,ending-start);
   //fclose(f1);

     ofstream myfile;
     myfile.open (resultFile,std::ios_base::app);
     myfile <<setprecision(30)<<treeId<<" "<<appResult<<" "<<overCounting<<" "<<ending-start<<endl;
     myfile.close();

   //cout<<sizeof(double)<<endl;
   //computeOvercountingFactor();
   //cout<<sizeof(bool)<<endl;
   //for (i=0;i<300;i++)
   //cout<<"num "<<setprecision(21)<<c.numColorings[c.num-1]<<" "<<ran()<<" "<<RANDMAX<<endl;


    deleteGraph(&g,n);
   // cout<<"error 1"<<endl;
    deleteTree(&t,m);
   // cout<<"error 2"<<endl;
    deleteColorings(&c,n,m);
   // cout<<"error 3"<<endl;
    for (i=0;i<=m;i++)
        free(subset[i]);
    //cout<<"error 4"<<endl;
    free(subset);
    //cout<<"error 5"<<endl;
    free(samplingTree);
    //cout<<"error 6"<<endl;
    free(samplingLoc);
    //cout<<"error 7"<<endl;
    for (i=0;i<n;i++)
        for (j=0;j<m;j++)
            free(num[i][j]);
    for (i=0;i<n;i++)
        free(num[i]);
    free(num);

    free(map);
    free(marked);

    //for (i=0;i<n;i++)
    //    free(edge[i]);
    //free(edge);

    //cout<<"error here";
    //system("PAUSE");
    return EXIT_SUCCESS;
}
