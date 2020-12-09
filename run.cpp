#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include <complex>
#include <fstream>
#include <list>
#include <algorithm>
#include <stack>
#include "YLM.h"
#include "read.h"


using namespace std; 
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//  unsigned seed = 1997;
std::mt19937_64 prng(seed); 
std::uniform_real_distribution<double> udist(0.0,1.0);

const double r_cut = 2.5; 
double lj_pair(double r_ij); 
double local_energy(int i);  
double total_energy();  
bool mcmove(); 
bool mcexc();
void init_pos();  
double dist_PBC(int i1, int i2); 
vector<double> dist_dir(int i1, int i2);
void calculate_complexQLM_6();
void if_solid();
void dot_q(int i, int j);

//  double T = 0.45;
//  double beta = 1/T;
//  int N = 50; 
//  double L = 44; 

double T = 0.45;
double beta = 1/T;
int N = 95; 
//double rho = 0.0025; 
//double L = pow(N/rho, 0.3333333); 
double L = 44.0;
int nd = 3; 
double V = pow(L, nd);
//double rho = N/V;
double P = 0.0005;
double mu = log(P/T)*T;
vector<double> pos(nd*N);
vector<vector<double>> list_q(N, vector<double>(26));
double total_e; 
vector<int> cluster_label; 

vector<int> particle_label;
vector<vector<int>> cluster_set; 
vector<vector<int>> allcluster_set;

void gcmc(); 
int find_cluster();
int find_allcluster();
/*----------------------------------------*/ 
int main ( int argc, char * argv[] ) {
  ofstream out_N("OP.out");
  ofstream check_pos("pos.out");

//    vector<double> pos;
  pos = getData("initial600.txt");
  int Npart = pos.size()/nd ;
  N = Npart; 

  int n = 3; // n is the n in Fig.2 of J. Phys.: Condens. Matter 21 (2009) 463102 
  int N0 = 800; 
  int n0 = 0; 
  int dt = 1; 
  vector<vector<double> > config_0; 
//    prng.seed(1213); 
  
  int n_measure = 0; 
  int n_step = 0; 
  int lambda0 = 4; 
  while(n0 < N0)
  {
    gcmc(); 
    n_step++;
    if(n_step%1000==0) // sample every 100 steps
    {
	    int lambda = find_allcluster();
	    cerr << n_step << "\t" << N << "\t" << lambda << "\t" << 0 << "\t" << 0 << endl; 
      if(lambda == lambda0) // this should be sampled with interval 
      {
	      n0++; 
        config_0.push_back(pos); 
      }
      n_measure++;
    } 
  }
  cout << "------------done with N0------------" << endl; 

  ofstream out_data("flux_P5e-4.out");
  double flux_0 =(double)n0/(double)n_measure;
  double flux_1; 
  cout << "end of interface 0, with flux0 = " << flux_0 << endl;
  cout << "number of stored config = " << config_0.size() << endl; 
  out_data << lambda0 << "\t" << flux_0 << endl;
  vector<vector<double>> config_i = config_0;
  vector<vector<double>> config_j;
  dt = 1;    
  int d_lambda = 4; 
  int n_itf = 1;
  vector<double> fluxes(n_itf); //itf indexes interface
  int lambdaf = lambda0;   
  for(int itf = 0; itf < n_itf; itf++)
  {
    std::uniform_int_distribution<int> uintdist(0, config_i.size()-1);
    lambdaf += d_lambda;
    int n_succ = 0; 
    int n_fail = 0; 
    int tot_trial = 0; 
    config_j = vector<vector<double> >(); 
    while(n_succ < 20)
    {
      int n = uintdist(prng); //??? no drawing new config here???
      pos = config_i[n]; // add: drawing from the collections 
      N = pos.size()/nd; 
      int n_step = 0;
      //int lambda = -1;
      int lambda = find_allcluster();
      while(((lambda > lambda0) && (lambda < lambdaf)) || (n_step == 0))
      {
        gcmc(); 

        n_step++; 
        lambda = find_allcluster();
//  	      if(n_step % 1 == 0)
//         	{
        cerr << n_step << "\t" << N << "\t" << lambda << "\t" << n_succ << "\t" << n_fail << endl; 
//  	      } 
      }
      tot_trial++; // sanity check: tot_trial != n_fail + n_succ
      if(lambda <= lambda0) // has to go back to the 0th interface
      {
        n_fail++; 
        cout << "FAIL :(" << endl; 
        cout << "lambda = " << lambda << ", n_fail: " << n_fail << endl; 
      }
      else if(lambda >= lambdaf)
      {
        n_succ++; 
        config_j.push_back(pos);
        cout << "SUCC :)" << endl; 
        cout << "lambda = " << lambda << ", n_succ: " << n_succ << endl; 
        ofstream succ_fs("succ_" + to_string(n_succ)); 
        for(int i = 0; i < pos.size(); i++)
        {
          succ_fs << pos[i] << endl;   
        }
      }
      else 
      {
        cout << "BUG IN LAMBDA1" << endl; 
      }
    }
    config_i = config_j;
    flux_1 = double(n_succ)/double(n_succ+n_fail); 
    fluxes.push_back(flux_1);
    cout << "n_succ: " << n_succ << endl; 
    cout << "n_fail: " << n_fail << endl; 
    cout << "flux_1: " << flux_1 << endl;
 	/* print out the flux */
    out_data << lambdaf << "\t" << n_succ << "\t" <<  n_fail << "\t" << tot_trial << "\t" << flux_1 << endl;    
  }
 

 return 0;
}
/*----------------------------------------*/ 
void gcmc()
{
  int swap_rate = 20; 
  std::uniform_int_distribution<int> uintdist(0, swap_rate);
  int randint = uintdist(prng);
  if(randint<swap_rate){
    mcmove();
  }
  else{
    mcexc();
  }
}
/*----------------------------------------*/ 
double lj_pair(double r_ij)
{
  double u_ij = 0; 
  if(r_ij > r_cut){
    u_ij = 0; 
  }
  else{
    double inv_rij = 1.0/r_ij;
    double inv_rij_2 = inv_rij * inv_rij;
    double inv_rij_6 = inv_rij_2 * inv_rij_2 * inv_rij_2;
    double inv_rij_12 = inv_rij_6 * inv_rij_6;
    u_ij   =  4.0 * (inv_rij_12 - inv_rij_6);   
  }

  return u_ij; 
}
/*----------------------------------------*/ 
double local_energy(int i) 
{
  int nPart = pos.size()/nd;
  double energy = 0.0; 
  for(int j = 0; j < nPart; j++)
  {
    if(j != i)
    {
      double dr = dist_PBC(i, j); 
      energy += lj_pair(dr); 
    }
  }

  return energy; 
}
/*----------------------------------------*/ 
double total_energy() 
{
  double energy = 0.0; 
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < i; j++)
    {
      double dr = dist_PBC(i, j); 
      energy += lj_pair(dr); 
    }
  }

  return energy; 
}
/*----------------------------------------*/ 
void init_pos() 
{
  int i,ix,iy,iz;
  
  int n3=2;
  /* Find the lowest perfect cube, n3, greater than or equal to the
     number of particles */
  while ((n3*n3*n3)<N) n3++;

  ix=iy=iz=0;
  /* Assign particle positions */
  for (i=0;i<N;i++) 
  {
    pos[i*nd+0] = ((double)ix+0.5)*L/n3; 
    pos[i*nd+1] = ((double)iy+0.5)*L/n3; 
    pos[i*nd+2] = ((double)iz+0.5)*L/n3; 
    ix++;
    if (ix==n3) 
    {
      ix=0;
      iy++;
      if (iy==n3) 
      {
        iy=0;
        iz++;
      }
    }
  }

  //cerr << "cluster in init: " << find_cluster() << endl; 
  //cerr << "n3: " << n3 << endl; 
  return; 
}
/*----------------------------------------*/ 
double dist_PBC(int i1, int i2)
//compute the distance between i1 and i2 in PBC 
{
  double d = 0.0;
  for (int k = 0; k < nd; k++)
  {
    // dr is in reduced space [-0.5 , 0.5 ) ^3
    int a_1 = i1*nd+k;
    int a_2 = i2*nd+k; 
    double dr = (pos[a_1] - pos[a_2])/L;
    dr = dr - rint(dr); // wrap back to MIC
    dr *= L;    // map to real space
    d += dr * dr; 
  }
  d = sqrt(d); 

  return d; 
}
/*-------------------------------------*/
bool mcmove()
{ 
  std::uniform_int_distribution<int> uintdist(0, N-1);
  int index = uintdist(prng);
  //Ccout << "index: " << index << endl; 
  
  vector<double> pos_old(nd);
  double eno = local_energy(index);
  for(int k=0; k<nd; k++)
  {
    double r = 0.4*(udist(prng)-0.5);
    int a_1 = index*nd + k;
    pos_old[k] = pos[a_1];
    pos[a_1] = fmod(pos[a_1]+r, L);
  }

  double enn = local_energy(index);
  bool acpt = false;
  if(udist(prng) < exp(-beta*(enn-eno))){
    acpt = true;
    total_e += enn-eno; 
  }
  else{
    for(int k=0; k< nd;k++){
      int a_2 = index*nd + k;
      pos[a_2] = pos_old[k];
    }
  }
  return acpt;
}
/*---------------------*/ 
bool mcexc()
{
  std::uniform_int_distribution<int> uintdist(0, N-1);

  double zz = exp(beta* mu);
  bool acpt = false;
  if(udist(prng)<0.5)
  {
    // try to delete
    if(N==0){
      return false;
    }
    else{
      int index = uintdist(prng);
      double enno = local_energy(index);
      double lnarg = log(N)-log(zz*L*L*L) + beta*enno;  
      double lnrand = log(udist(prng));
      if(lnrand < lnarg){
        pos.erase(pos.begin()+index*nd, pos.begin()+index*nd+nd);
        N = N-1;
        acpt = true;
	//cout << "delete" << endl;
      }
    }
  }
  else
  {
  // insert
    vector<double> add_pos(nd);
    for(int k=0;k<nd;k++){
      add_pos[k] = L* udist(prng);
      pos.push_back(add_pos[k]);
    }
    // add the new particle to the end of pos
    double enadd = local_energy(N);
   // cout << "del E for insertion: " << enadd << " Npart = "<< pos.size()/nd << "  N = " << N << endl;
   // cout << "position: " << add_pos[0] << "\t" << add_pos[1] << "\t" << add_pos[2] << endl;
    double lnarg = log(zz*L*L*L) - log(N+1) -beta*enadd;
    double lnrand = log(udist(prng));
    if (lnrand < lnarg){
      N = N+1;
      acpt =true;
    }
    else{
      pos.erase(pos.begin()+N*nd, pos.begin()+N*nd+nd);
    }
  }     
       
  return acpt; 
}

/*-----------------------------*/
bool if_neighbor(int id_node, int id_check)
{
  bool ifnb = false;
  if(dist_PBC(id_node, id_check)<1.5){
    ifnb = true;
  }
  return ifnb;
}
/*----------------------------------------*/
double dist_PBC(vector<double>& pos, int i1, int i2)
//compute the distance between i1 and i2 in PBC
{
  double d = 0.0;
  for (int k = 0; k < nd; k++)
  {
    int a_1 = i1*nd+k;
    int a_2 = i2*nd+k;
    double hL = L/2.0;
    double dr = (pos[a_1] - pos[a_2]);
    if(dr>hL){
      dr = dr - L;
    }
    if(dr<-hL){
      dr = dr + L;
    }

    d += dr * dr;
  }
  d = sqrt(d);

  return d;
}
/*-----------------------------*/
bool if_neighbor(vector<double>& pos, int id_node, int id_check){
  bool ifnb = false;
  if(id_node!= id_check && dist_PBC(pos, id_node, id_check)<1.5){
    ifnb = true;
  }
  return ifnb;
}
/*----------------------------rewrite the cluster OP*/
vector<int> n_neighbors(vector<double>& pos, int index)
{
  vector<int> nb_list;
  int Npart = pos.size()/nd;
  for(int j=0; j < Npart; j++)
  {
    if(dist_PBC(pos, index, j)<1.5 && j!=index)
    {
      nb_list.push_back(j);
    }
  }
  return nb_list;
}
/*--------------------------*/
vector<double> dist_dir(vector<double>& pos, int i1, int i2){
  //displacement, not normalized
  vector<double> dvect_r(nd);
  double hL = L/2.;
  for(int k=0; k <nd; k++)
  {
	  int a1 = i1*nd+k;
	  int a2 = i2*nd+k;
  	dvect_r[k] = pos[a1] - pos[a2];
    if (dvect_r[k] > hL){
      dvect_r[k] = dvect_r[k] - L;
    }
    if (dvect_r[k] < -hL){
      dvect_r[k] = dvect_r[k] + L;
    }
  }
  return dvect_r;
}
/*--------------------------*/
void calculate_complexQLM_6(vector<double>& pos, vector<vector<double>>& list_q){
  //nn = number of neighbors
  //list_q contains q_6m for each particle; q[0-12] = realq, q[13-25]=imagq;
  double realti,imgti;
  double realYLM,imgYLM;
  int N = list_q.size();
  for (int t= 0; t < N; t++)
  {
    // cout << t << endl;
    vector<int> nn_list = n_neighbors(pos, t);
    int nn = nn_list.size();
	  for (int mi = -6; mi < 7; mi++)
    {
      realti = 0.0;
      imgti = 0.0;
      for (int c = 0; c < nn; c++)
      {
  	    vector<double> dvec_r = dist_dir(pos, t, nn_list[c]);

        double dx = dvec_r[0];
  	    double dy = dvec_r[1];
        double dz = dvec_r[2];
  	    double theta;
        double phi;
      	double r;

        convert_to_spherical_coordinates(dx, dy, dz, r, phi, theta);

      	QLM(6, mi, theta, phi,realYLM, imgYLM);

		    realti += realYLM;
		    imgti += imgYLM;
	    }
//        realti = nn == 0 ? 0: realti;
//        imgti = nn == 0 ? 0: imgti;
      realti = nn == 0 ? 0 : realti/(double(nn));
      imgti = nn == 0? 0 : imgti/(double(nn)); //FIXME: don't know whether to divide nn.
      list_q[t][mi+6] = realti;
      list_q[t][13+mi+6] = imgti;
    }
  }
}
/*--------------------------*/
complex<double> dot_q(vector<vector<double>>& list_q, int i, int j){
  complex<double> product = 0;
  for(int mi = -6; mi < 7; mi++)
  {
    // cout << "list_q[i][mi+6]: " << list_q[i][mi+6] << endl;
    // cout << "list...: " << list_q[i][12+mi+6] << endl;
    // cout << "666: " << (list_q[i][mi+6], list_q[i][12+mi+6]) << endl;
    complex<double> q6i = complex<double>(list_q[i][mi+6], list_q[i][13+mi+6]);
    complex<double> q6j = complex<double>(list_q[j][mi+6], list_q[j][13+mi+6]);

    product = product + q6i * conj(q6j);
  }
  return product;
}

/*--------------------------*/
// q_lm is now normalized when calculating the second order parameter
/*--------------------------*/
// void normalized_complexQLM_6(vector<vector<double>>& list_q, int& N){
//
//   for(int i =0; i < N; i++){
//     double norm_i = norm_q(list_q[i]);
//     for (int mi = -6; mi < 7; mi++){
//
//       list_q[i][mi+6] =  norm_i == 0 ? 0 : list_q[i][mi+6]/ norm_i;
//       list_q[i][13+mi+6] = norm_i == 0 ? 0 :list_q[i][13+mi+6]/ norm_i;
//
//     }
//   }
// }
/*--------------------------*/
double norm_q(vector<double> vec_q) {
    double normq = 0.;
    int dim_q = vec_q.size()/2; // double the size bc of real and imag part
    if(dim_q != 13){
        cout << "BUG from q6" << endl;
    }
    for (int mi = -6; mi < 7; mi++){

      normq = normq + vec_q[mi+6] * vec_q[mi+6] + vec_q[mi+19] * vec_q[mi+19];
    }
    normq = sqrt(normq);
  return normq;
}
///*--------------------------*/
// Q6 is not used in OP; just for sanity check
// void compute_Q6(vector<vector<double>>& list_q, int& N){
// // construct a rotational invariant Q_l out of average of q_lm over all particles
//     // calculate the average Q_lm = \sum(q_lm(i)) /N
//     vector<double> avgQ_lm(26);
//     // cout << N << endl;
//     for(int mi = -6; mi < 7; mi++){
//         for (int i=0; i < N; i++){
//              // cout << list_q[i][mi+6] << endl;
//              avgQ_lm[mi+6] = avgQ_lm[mi+6] + list_q[i][mi+6]/N;
// 	           avgQ_lm[13+mi+6] = avgQ_lm[13+mi+6] + list_q[i][13+mi+6]/N;
//        }
//     }
//     // calculate the Q_l = sqrt(4 pi/(2l+1)) \sum_m |Q_lm|^2
//
//     double Q_6 = sqrt(4*PI/13.0) * norm_q(avgQ_lm);
//     // cout << "norm(avg_q) = " << norm_q(avgQ_lm) << endl;
//     cout << "Q_6 = " << Q_6 << endl;
// }
// //*--  ------------------------*/
// // 1. compute q_lm for all particles
// // 2. find all solid-like particles
///////  solid-like if having more than 5 bonds
// // 3. find all solidlike cluster and return the size of the largest cluster
// //*--------------------------*/
/*----------------------------the second cluster OP*/
int find_solidcluster(vector<double>& pos, vector<vector<double>>& list_q){
  // particle_label = vector<int>(N, -1);
  // allcluster_set = vector<vector<int>>(0);
  int N = list_q.size();
  int cnt_sp = 0;
  vector<int> solid;

  for(int i=0; i<N; i++){
    int n_bond = 0;
    vector<int> nl = n_neighbors(pos, i);
    for(int k=0; k<nl.size(); k++){
      int j = nl[k];
      if (j!=i) {
        complex<double> product = dot_q(list_q, i, j);
        complex<double> ni = dot_q(list_q, i, i);
        complex<double> nj = dot_q(list_q, j, j);
        // double prod = real(product);
        double prod = real(product)/sqrt(real(ni) * real(nj));

        if(prod > .65){
          n_bond++;
        }
      }
    }
    // cout << "particle " << i << "with # bond = " << n_bond << endl;
    if (n_bond> 5){
      solid.push_back(i);
      cnt_sp++;
      }
    }
    // cout << "found num solidlike = " << cnt_sp << endl;

  // do cluster search among all the solid particles
  vector<int> particle_label(cnt_sp, -1);

  int max_cluster = -1;
  int cluster_label = 0;

  for(int start=0; start < cnt_sp; start++)
  {
    if(particle_label[start] == -1)
    {
      cluster_label++;
      int cluster_size = 0;
      vector<bool> hasnot_entered(cnt_sp, true);
      vector<int> new_cluster;
      stack<int> stack;
      stack.push(start);
      new_cluster.push_back(start);
      particle_label[start] = cluster_label;
      hasnot_entered[start] = false;

      while(!stack.empty())
      {
        int origin = stack.top();
        stack.pop();

        for(int j = 0; j < cnt_sp;j++)
        {
          if( hasnot_entered[j] &&  if_neighbor(pos, solid[origin], solid[j]))
//            if(particle_label[j] != -1 &&  if_neighbor(pos, solid[origin], solid[j]))
          {
            new_cluster.push_back(j);
            particle_label[j] = cluster_label;
            stack.push(j);
	          hasnot_entered[j] = false;
          }
        }
      }


      // cout << "cluster label: " << cluster_label << endl;
      // cout << " with member: " << endl;
      // for(int i=0; i< cluster_size; i++){
      //   cout << new_cluster[i] << endl;
      // }

      cluster_size = int(new_cluster.size());
      max_cluster = max(max_cluster, int(new_cluster.size()));

      // allcluster_set.push_back(new_cluster);

      //cout << "largest cluster_size ="<< cluster_size << endl;
      //cout << "begin w/ particle:" << start << "found size="<< cluster_size << " with label" << cluster_label << endl;
    }
  }
//    cout << "cnt_sp: " << cnt_sp << endl; 
  return max_cluster;
}

/*--------------------------*/
int find_allcluster(){

  vector<vector<double>> list_q(N, vector<double>(26));
  calculate_complexQLM_6(pos,list_q);
  int sec_op = find_solidcluster(pos, list_q);

  return sec_op; 
}
/*--------------------------*/
vector<int> n_neighbors(int index)
{
  vector<int> nn;
  for(int j=0;j<N;j++)
  {
    if(dist_PBC(index,j)<1.5 && j!=index)
    {
      nn.push_back(j);
    }
  }
  return nn;
}
/*--------------------------*/
vector<double> dist_dir(int i1, int i2){
  //displacement, not normalized
  vector<double> dvect_r(nd);
  double hL = L/2.;
  for(int k=0; k <nd; k++)
  {
	  int a1 = i1*nd+k;
	  int a2 = i2*nd+k;
  	dvect_r[k] = pos[a1] - pos[a2];
    if (dvect_r[k] > hL){
      dvect_r[k] = dvect_r[k] - L;
    }
    if (dvect_r[k] < -hL){
      dvect_r[k] = dvect_r[k] + L;
    }
  }
  return dvect_r;
}
/*--------------------------*/
void calculate_complexQLM_6()
{
  //nn = number of neighbors
  //turn realq/imagq to global variable?
  //vector<vector<double>> list_q(N); 
  //list_q (global variable) contains q6 for each particle; q[0-12] = realq, q[13-25]=imagq;
  int nn;
  double realti,imgti;
  double realYLM,imgYLM;
  // nop = parameter.nop;
  for (int t= 0; t<N; t++)
  {
    vector<int> nn_list = n_neighbors(t);
    int nn = nn_list.size();
	  for (int mi = -6; mi < 7; mi++)
    {
      realti = 0.0;
      imgti = 0.0;
      for (int c = 0; c < nn; c++)
      {
  	    vector<double> dvec_r = dist_dir(t, nn_list[c]);
        double dx = dvec_r[0];
  	    double dy = dvec_r[1];
        double dz = dvec_r[2];
  	    double theta;
        double phi;
      	double r;
        convert_to_spherical_coordinates(dx, dy, dz, r , phi, theta);
        if(mi == -6)
        {
//            cout << "t: " << t << endl; 
//            cout << "nn_list[c]: " << c << ", " << nn_list[c] << endl; 
//            cout << "dx: " << dx << endl; 
//            cout << "dy: " << dy << endl; 
//            cout << "dz: " << dz << endl; 
//            cout << "theta: " << theta << endl; 
//            cout << "phi: " << phi << endl; 
        }

      	QLM(6, mi, theta, phi,realYLM, imgYLM);
		    realti += realYLM;
		    imgti += imgYLM;
	    }

      //cout << "t, m:" << t << ", " << mi << endl; 
      //cout << "realit: " << realti << endl; 
      //cout << "imagti: " << imgti << endl; 
      
      realti = nn == 0 ? 0 : realti/(double(nn));
      imgti = nn == 0? 0 : imgti/(double(nn));
      list_q[t][mi+6] = realti;
      list_q[t][12+mi+6] = imgti;
    }

//  atoms[ti].realq[4][mi+6] = realti;
//  atoms[ti].imgq[4][mi+6] = imgti;
//  no idea why realq[4], but mi+6 starts with 0}
  }
}
/*--------------------------*/
// complex<double> dot_q(int i, int j){
void dot_q(int i, int j){
  complex<double> product;

  for(int mi = -6; mi < 7; mi++)
  {
    cout << "list_q[i][mi+6]: " << list_q[i][mi+6] << endl; 
    cout << "list...: " << list_q[i][12+mi+6] << endl; 
    cout << "666: " << (list_q[i][mi+6], list_q[i][12+mi+6]) << endl; 
    complex<double> q6i = complex<double>(list_q[i][mi+6], list_q[i][12+mi+6]);
    complex<double> q6j = complex<double>(list_q[j][mi+6], list_q[j][12+mi+6]);
    product = product + q6i * conj(q6j);
    cout << real(q6i) << ", " << real(q6j) << endl;    
  }
  cout << "compute dot product: "<< real(product) << "," << imag(product) << endl;
  cout << "dist: " << dist_PBC(i, j) << endl; 

}
/*--------------------------*/
void if_solid(int i, int j){
	// q6(i) is a (1*13) vector
	// compute q6(i) dot q6(j)

}
