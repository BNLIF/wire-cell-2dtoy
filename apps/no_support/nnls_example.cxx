
#include <lsp/nnls.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <sstream>

int main(){
  using namespace boost::numeric::ublas;
  int m = 6, n=5;
  matrix< double > A(n,m);
  vector< double > b(n);
  A(0,3) = 1; A(0,4) = 1; A(0,5) = 1;
  A(1,0) = 1; A(1,1) = 1; A(1,2) = 1;
  A(2,2) = 1; A(2,5) = 1;
  A(3,1) = 1; A(3,4) = 1;
  A(4,0) = 1; A(4,3) = 1;
  
  b(0) = 28.7355;
  b(1) = 7.86395;
  b(2) = 18.0824;
  b(3) = 15.8977;
  b(4) = 2.61937;

  lsp::nnls< matrix< double >, vector< double >  > nnls( A, b );
  vector< double > x(m);
  matrix< double > cov(m,m);
  nnls.solve( x , cov );

  for (int i=0;i!=m;i++){
    std::cout << x(i) << std::endl;
  }
  
  vector<double> r = prod(A,x) -b;
  std::cout << "Chi2: " <<  inner_prod(r,r) << std::endl;
}
