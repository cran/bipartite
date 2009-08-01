// key point is lines 147 ff.: algorithm to find min matrix


#include <iostream>
#include <iterator>
#include <vector>
#include <map>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <cmath>

#include "crossplatform.hpp"


#include "matrix.hpp"
//#include "rxc_enumwrapper.hpp"


using namespace std;
using namespace TNT;

double minimum=1000000.0;
double maximum=0.0;

int sum(const vector<int> &n) {
  int tmp=0;
  for (vector<int>::const_iterator iter=n.begin(); iter!=n.end(); iter++) {
    tmp+=*iter;
  }
  return tmp;
}


int getMaxIndex(const vector<int> &x) {
  int max=0;
  int index=-1;
  for (uint i=0; i<x.size(); ++i) {
    if (x[i]>max) { index=i; max=x[i]; }
  }

  return index;
}

double getS( const Array2D<int> &a ) 
{
  double stat=0;
  vector<int> n,m;
  getRowColSums( n, m, a );
  int total=sum(n);
  for (int i=0; i<a.dim1(); ++i) {
    for (int j=0; j<a.dim2(); ++j) {
      if (a[i][j]>0) {
        stat+=(double)a[i][j]/(double)total*log((double)a[i][j]/(double)total);
      }
    }
  }
  return -stat;
}

double getT( const Array2D<int> &a ) 
{
  double stat=0;
  for (int i=0; i<a.dim1(); ++i) {
    for (int j=0; j<a.dim2(); ++j) {
      if (a[i][j]>0) {
        stat+=(double)a[i][j]*log((double)a[i][j]);
      }
    }
  }
  return stat;
}

void getRowColSums( vector<int> &n, vector <int> &m, const Array2D<int> &a )
{

  n.clear();
  m.clear();
  n.resize(a.dim1(),0);
  m.resize(a.dim2(),0);
  for (int i=0; i<a.dim1(); ++i) {
    for (int j=0; j<a.dim2(); ++j) {
      n[i]+=a[i][j];
      m[j]+=a[i][j];
    }
  }
  
}

std::ostream &operator<<(std::ostream &ostr, vector<int> v) {
  copy(v.begin(), v.end(), ostream_iterator<int>(ostr,"\t"));
  return ostr;
}

Array2D<int> getMaxArray(const vector<int> &n, const vector<int> &m) {
  vector<int> ntmp=n, mtmp=m;
  Array2D<int> workingarray( n.size(), m.size(), 0  );
  // Fill all maximally available matrix elements
  int indmaxn=getMaxIndex(ntmp);
  int indmaxm=getMaxIndex(mtmp);
  while ((indmaxm>=0)&&(indmaxn>=0)) {
    // Find indmaxn and indmaxm for indexes of maximal border values
    assert(((indmaxm<0)&&(indmaxn<0))
           ||((indmaxm>=0)&&(indmaxn>=0)));
 
    // Take the minimum as the new entry in the matrix, and adjust the
    // remaining row and colsums that have to be distributed along the
    // array.
    uint entry=std::min(ntmp[indmaxn],mtmp[indmaxm]);
    if (entry!=0) {
      workingarray[indmaxn][indmaxm]=entry;
      ntmp[indmaxn]-=entry;
      mtmp[indmaxm]-=entry;
    }
    indmaxn=getMaxIndex(ntmp);
    indmaxm=getMaxIndex(mtmp);
  } 
  return workingarray;
}


Array2D<int> getMinArray(const vector<int> &n, const vector<int> &m) {
  int N=0;
  for (uint i=0; i<n.size(); ++i) {
    N+=n[i];
  }
  vector<int> ntmp, mtmp;
  Array2D<int> workingarray( n.size(), m.size(), 0 );
  multimap<int,pair<uint,uint> > entries;
  for (uint i=0; i<n.size(); ++i) {
    for (uint j=0; j<m.size(); ++j) {
      workingarray[i][j]=(int)floor((double)n[i]*(double)m[j]/(double)N);
    }
  }

  getRowColSums(ntmp,mtmp,workingarray);

  bool finished;
  
  do {

    finished=1;
    for (uint i=0; i<ntmp.size(); ++i) {
      if (n[i]>ntmp[i]) finished=0;
    }
    for (uint j=0; j<mtmp.size(); ++j) {
      if (m[j]>mtmp[j]) finished=0;
    }
    
    if (!finished) {
      entries.clear();
      for (uint i=0; i<n.size(); ++i) {
        for (uint j=0; j<m.size(); ++j) {
          entries.insert(make_pair((int)workingarray[i][j],make_pair(i,j)));
        }
      }
      
      for (multimap<int,pair<uint,uint> >::iterator iter=entries.begin();
           iter!=entries.end(); ++iter) {
        uint i=iter->second.first;
        uint j=iter->second.second;
        if ((n[i]>ntmp[i]) && (m[j]>mtmp[j])) {
          workingarray[i][j]++;
          ntmp[i]++;
          mtmp[j]++;
        }
      }
    } 
    
  } while (!finished);
  

  return workingarray;
}


double getMax(const vector<int> &n, const vector<int> &m) {
  return getS(getMaxArray(n,m));
}

double getMin(const vector<int> &n, const vector<int> &m) 
{
  return getS(getMinArray(n,m));
}

//double getMin(const vector<int> &n, const vector<int> &m) 
//{
//  int N;
//  for (uint i=0; i<n.size(); ++i) {
//    N+=n[i];
//  }
//  vector<int> ntmp=n, mtmp=m;
//  Array2D<int> workingarray( n.size(), m.size(), 0 );
//  
//  int counter=N;
//  do {
//    for (uint i=0; i<ntmp.size(); ++i) {
//      for (uint j=0; j<mtmp.size(); ++j) {
//        if ((ntmp[i]>0) && mtmp[j]>0) {
//          workingarray[i][j]++;
//          ntmp[i]--;
//          mtmp[i]--;
//          counter--;
//          if (counter<=0) 
//            return getT(workingarray);
//        }
//      }
//    }
//  } while (counter>0);
//  
//  return getT(workingarray);
//}

// External declaration for the integrator

//void test( const Array2D<int> &a ) {
//  
//  vector<int> n,m;
//  getRowColSums(n,m,a);
//
//  double maxT=getMax(n,m);
//  double minT=getMin(n,m);
//
//  // This ugly code is to interface the fortran library.
//  int r,c,*nf,*mf,ifault;
//  r=n.size();
//  c=m.size();
//  nf=&n[0];
//  mf=&m[0];
//
//  enumber_( &r, &c, 
//       nf, 
//       mf, 
//       &ifault );    
//  
//  if (ifault) {
//    cout << a << "IFAULT:" << ifault << endl;
//  }
//
//  
//  cout << minimum << "\t" << maximum;
//  if (maximum!=maxT) {
//    cout << a << "maximum:" << maxT << "!=" << maximum << endl;
//  }
//  if (minimum!=minT) {
//    cout << a << "minimum:" << minT << "!=" << minimum << endl;
//  }
//}


vector<double> KLD(vector<int> n, vector<int> m, 
                   Array2D<int> a, int total, 
                   bool transpose)
{
  vector<double> tmp(n.size(),0.0);
  for (uint i=0; i<n.size(); ++i) {
    // KD
    for (uint j=0; j<m.size(); ++j) {
      double fij=(double)(transpose?a[j][i]:a[i][j])/(double)n[i];
      double gj=(double)m[j]/(double)total;
      if (fij>0.0) {
        tmp[i]+=fij*log(fij/gj);
      }
    }
  }
  return tmp;
}


struct index_count { 
  index_count(uint i, int c) : index(i), count(c) {}
  uint index;
  int count;
  bool operator<(const index_count& a) const {
    return count<a.count;
  }
  bool operator<(const int& a) const {
    return count<a;
  }
};

class NotSmallerIndexCount {
public:
  NotSmallerIndexCount(int count) : 
    count_(count) {};

  bool operator()(const index_count &ic) const {
    return !(ic<count_);
  }

  int count_;
};


vector<double> maxKLD(vector<int> n, vector<int> m, 
                   int total)
{
  vector<double> maxkld(n.size(),0.0);

  for (uint i=0; i<n.size(); ++i) {
    maxkld[i]=log((double)total/(double)n[i]);
  }
  return maxkld;
}

bool operator<(const pair<double,uint>& a, const pair<double,uint>& b) {
    return a.first < b.first;
}

double kld(double p, double q) {
  if (p==0.0) return 0;
  return p*log(p/q);
}

void distribute_for_kld(uint sum, const int &n, const vector<int> &m, 
                        int &total, vector<int> &current ) {

  if (sum>=n) return;
  vector<pair<double,uint> > index;
  for (uint j=0; j<m.size(); ++j) {
    double delta=kld((double)(current[j]+1)/(double)n,((double)m[j]/(double)total))
      -kld((double)current[j]/(double)n,((double)m[j]/(double)total));
    index.push_back(pair<double,uint>(delta, j));
  }
  sort(index.begin(), index.end());
  current[index[0].second]++;
  distribute_for_kld(sum+1, n, m, total, current);
}

vector<double> minKLD(vector<int> n, vector<int> m, 
                   int total) 
{
  vector<double> minkld(n.size(),0.0);
  
  for (uint i=0; i<n.size(); ++i) {
    int sum=0;
    vector<int> current(m.size(),0);
    for (uint j=0; j<m.size(); ++j) {
      current[j]=(int)((double)n[i]*((double)m[j]/(double)total));
      sum+=current[j];
    }
    distribute_for_kld(sum, n[i], m, total, current );
    sum=0;
    for (uint j=0; j<m.size(); ++j) {
      sum+=current[j];
      minkld[i]+=kld((double)current[j]/(double)n[i],
                     ((double)m[j]/(double)total));
    }
    assert(sum==n[i]);
    assert(minkld[i]>=0.0);
  }
  return minkld;
}
