#include <Rcpp.h>
using namespace Rcpp;

double abss (double x){
  if (x<0) return -x; 
  else return x;
}

double max(double a, double b)
{
  if (a<b) 
    return b; 
  else 
    return a;
}


NumericVector dist(NumericVector x, double z, int n){
  NumericVector d(n);
  for(int i = 0; i < n; i++){
    d[i] = abss(x[i]-z);
  }
  return d;
}

int comp(double x, double y){
  if (x <= y) 
    return 1; 
  else 
    return 0;
}


int cal_aij(NumericVector x, NumericVector y, double x0, double y0, int n, int c)
{ 
  int m = 0;
  for(int i = 0; i < n; i++)
  {
    int temp=0;
    switch(c){
    case 0: temp = comp(x[i], x0)*comp(y[i], y0); break; // P11 Q11
    case 1: temp = comp(x[i], x0)*(1 - comp(y[i], y0)); break; // P12 Q12
    case 2: temp = (1 - comp(x[i], x0))*comp(y[i], y0); break; // P21 Q21
    case 3: temp = (1 - comp(x[i], x0))*(1 - comp(y[i], y0)); break; // P22 Q22
    }
    m = m + temp;
  }
  return m;
}

// cal Pij and Qij
NumericVector cal_pq(NumericVector x, NumericVector y, int n, double z1, double z2, double da, double db){
  NumericVector res(4);
  // 先计算x和y到该点的距离
  NumericVector resa = dist(x, z1, n);
  NumericVector resb = dist(y, z2, n);
  
  //Rcout << "da: " << da << " db: " << db << "\n";
  //Rcout << "resa: " << resa[0] << resa[1] << resa[2] << resa[3] << "\n"; 
  //Rcout << "resb: " << resb[0] << resb[1] << resb[2] << resb[3] << "\n";
  
  res[0] = cal_aij(resa, resb, da, db, n, 0);  // P11 Q11
  res[1] = cal_aij(resa, resb, da, db, n, 1);  // P12 Q12
  res[2] = cal_aij(resa, resb, da, db, n, 2);  // P21 Q21
  res[3] = cal_aij(resa, resb, da, db, n, 3);  // P22 Q22
  
  return res;
}



// 计算T2
double cal_t (NumericVector p, NumericVector q, int nx, int ny)
{
  double t = 0;
  double r1;
  double r2;
  for(int i=0; i<4; i++)
  {
    r1 = (p[i] + q[i])*nx/(nx+ny);
    r2 = (p[i] + q[i])*ny/(nx+ny);
    if(r1 != 0){
      if(r2 != 0){
        t = t + (p[i] - r1)*(p[i] - r1) / r1 + (q[i] - r2)*(q[i] - r2) / r2;
      }
    }
    
  }
  return t;
}


// 计算 MAC
double macc2(NumericVector xi1, NumericVector xi2, NumericVector yj1, NumericVector yj2, NumericVector z1, NumericVector z2, int nx, int ny, int k){
  double mac = 0;
  
  int dn = k;               // 用于计算距离的点的个数
  
  //Rcout << dn << "\n";
  
  for(int i = 0; i < dn; i++){
    //for(int j = 0; j < dn; j++){
    for(int j = i+1; j < dn; j++){  
      double da = abss(z1[i] - z1[j]);  // Axiyj
      double db = abss(z2[i] - z2[j]);  // Bxiyj
      
      //Rcout << "i=" << i << " j=" << j << ":" << " da=" << da << " db=" << db << "\n";
      
      NumericVector p(4);           // Pij
      NumericVector q(4);           // Qij
      
      p = cal_pq(xi1, xi2, nx, z1[i], z2[i], da, db);
      
      //Rcout << "p: " << p[0] << p[1] << p[2] << p[3] << "\n";
      
      q = cal_pq(yj1, yj2, ny, z1[i], z2[i], da, db);
      
      //Rcout << "q: " << q[0] << q[1] << q[2] << q[3] << "\n";
      
      double t = cal_t(p, q, nx, ny);   // 计算t
      
      //Rcout << "t: " << t << "\n";
      mac = max(mac, t);
    }
  }
  
  
  return mac;
}



// [[Rcpp::export]]
double mac2(NumericMatrix X, NumericMatrix Y, NumericMatrix Z, int nx, int ny, int k)  //Main function
{
  NumericVector xi1(nx);  // X[0]   
  NumericVector xi2(nx);  // X[1]
  NumericVector yj1(ny);  // Y[0]
  NumericVector yj2(ny);  // Y[1]
  for(int i=0; i<nx; i++){
    xi1[i] = X(i,0);
    xi2[i] = X(i,1);
  }
  for(int i=0; i<ny; i++){
    yj1[i] = Y(i,0);
    yj2[i] = Y(i,1);
  }
  
  if(nx<k) k=nx;
  if(ny<k) k=ny;
  
  NumericVector z1(k);
  NumericVector z2(k);
  for(int i=0; i<k; i++){
    z1[i] = Z(i,0);
    z2[i] = Z(i,1);
  }
  
  // define new function to cal mac
  double res = macc2(xi1,xi2,yj1,yj2,z1,z2,nx,ny,k);
  
  return res;
}


