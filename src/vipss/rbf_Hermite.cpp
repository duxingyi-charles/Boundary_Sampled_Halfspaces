#include "rbfcore.h"
#include "utility.h"
#include "Solver.h"
#include <armadillo>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <queue>
#include "readers.h"
//#include "mymesh/UnionFind.h"
//#include "mymesh/tinyply.h"

typedef std::chrono::high_resolution_clock Clock;
double randomdouble() {return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);}
double randomdouble(double be,double ed) {return be + randomdouble()*(ed-be);	}

void RBF_Core::NormalRecification(double maxlen, vector<double>&nors){


    double maxlen_r = -1;
    auto p_vn = nors.data();
    int dim = point_dimension;
    int  np = nors.size()/dim;
    if(1){
        for(int i=0;i<np;++i){
            maxlen_r = max(maxlen_r,MyUtility::normVec(p_vn+i*dim,dim));
        }

        cout<<"maxlen_r: "<<maxlen_r<<endl;
        double ratio = maxlen / maxlen_r;
        for(auto &a:nors)a*=ratio;
    }else{
        for(int i=0;i<np;++i){
            MyUtility::normalize(p_vn+i*3);
        }

    }




}

bool RBF_Core::Write_Hermite_NormalPrediction(string fname, int mode){


//    vector<uchar>labelcolor(npt*4);
//    vector<uint>f2v;
//    uchar red[] = {255,0,0, 255};
//    uchar green[] = {0,255,0, 255};
//    uchar blue[] = {0,0,255, 255};
//    for(int i=0;i<labels.size();++i){
//        uchar *pcolor;
//        if(labels[i]==0)pcolor = green;
//        else if(labels[i]==-1)pcolor = blue;
//        else if(labels[i]==1)pcolor = red;
//        for(int j=0;j<4;++j)labelcolor[i*4+j] = pcolor[j];
//    }
    //fname += mp_RBF_METHOD[curMethod];

//    for(int i=0;i<npt;++i){
//        uchar *pcolor = green;
//        for(int j=0;j<4;++j)labelcolor[i*4+j] = pcolor[j];
//    }

    vector<double>nors;
    if(mode ==0)nors=initnormals;
    else if(mode == 1)nors=newnormals;
    else if(mode == 2)nors = initnormals_uninorm;
    NormalRecification(1.,nors);

    //for(int i=0;i<npt;++i)if(randomdouble()<0.5)MyUtility::negVec(nors.data()+i*3);
    //cout<<pts.size()<<' '<<f2v.size()<<' '<<nors.size()<<' '<<labelcolor.size()<<endl;
    //writePLYFile(fname,pts,f2v,nors,labelcolor);

//    writeObjFile_vn(fname,pts,nors);
    if (point_dimension == 3) {
        writePLYFile_VN(fname,pts,nors);
    } else { //2D
        writePLYFile_VN_2D(fname,pts,nors);
    }

    return 1;
}



//void RBF_Core::Set_HermiteRBF(vector<double>&pts){
//
//    cout<<"Set_HermiteRBF"<<endl;
//    //for(auto a:pts)cout<<a<<' ';cout<<endl;
//    isHermite = true;
//
//    a.set_size(npt*4);
//    M.set_size(npt*4,npt*4);
//    double *p_pts = pts.data();
//
//    //M00
//    for(int i=0;i<npt;++i){
//        for(int j=i;j<npt;++j){
//            M(i,j) = M(j,i) = Kernal_Function_2p(p_pts+i*3, p_pts+j*3);
//        }
//    }
//
//
//    //if(User_Lamnbda!=0)for(int i=0;i<npt;++i)M(i,i) += User_Lamnbda;
//
//    //M01
//    double G[3];
//    for(int i=0;i<npt;++i){
//        for(int j=0;j<npt;++j){
//
//            Kernal_Gradient_Function_2p(p_pts+i*3, p_pts+j*3, G);
//            //            int jind = j*3+npt;
//            //            for(int k=0;k<3;++k)M(i,jind+k) = -G[k];
//            //            for(int k=0;k<3;++k)M(jind+k,i) = G[k];
//
//            for(int k=0;k<3;++k)M(i,npt+j+k*npt) = G[k];
//            for(int k=0;k<3;++k)M(npt+j+k*npt,i) = G[k];
//
//        }
//    }
//
//    //M11
//    double H[9];
//    for(int i=0;i<npt;++i){
//        for(int j=i;j<npt;++j){
//
//            Kernal_Hessian_Function_2p(p_pts+i*3, p_pts+j*3, H);
//            //            int iind = i*3+npt;
//            //            int jind = j*3+npt;
//            //            for(int k=0;k<3;++k)
//            //                for(int l=0;l<3;++l)
//            //                    M(jind+l,iind+k) = M(iind+k,jind+l) = -H[k*3+l];
//
//            for(int k=0;k<3;++k)
//                for(int l=0;l<3;++l)
//                    M(npt+j+l*npt,npt+i+k*npt) = M(npt+i+k*npt,npt+j+l*npt) = -H[k*3+l];
//        }
//    }
//
//    //cout<<std::setprecision(5)<<std::fixed<<M<<endl;
//
//    bsize= 4;
//    N.zeros(npt*4,4);
//    b.set_size(4);
//
//    for(int i=0;i<npt;++i){
//        N(i,0) = 1;
//        for(int j=0;j<3;++j)N(i,j+1) = pts[i*3+j];
//    }
//    for(int i=0;i<npt;++i){
//        //        int ind = i*3+npt;
//        //        for(int j=0;j<3;++j)N(ind+j,j+1) = 1;
//
//        for(int j=0;j<3;++j)N(npt+i+j*npt,j+1) = -1;
//    }
//
//    //cout<<N<<endl;
//    //arma::vec eigval = eig_sym( M ) ;
//    //cout<<eigval.t()<<endl;
//
//
//    if(!isnewformula){
//        cout<<"start solve M: "<<endl;
//        auto t1 = Clock::now();
//        if(isinv)Minv = inv(M);
//        else {
//            arma::mat Eye;
//            Eye.eye(npt*4,npt*4);
//            Minv = solve(M,Eye);
//        }
//        cout<<"solved M: "<<(invM_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;
//
//        t1 = Clock::now();
//        if(isinv)bprey = inv_sympd(N.t() * Minv * N) * N.t() * Minv;
//        else {
//            arma::mat Eye2;
//            Eye2.eye(bsize,bsize);
//            bprey = solve(N.t() * Minv * N, Eye2) * N.t() * Minv;
//        }
//        cout<<"solved bprey "<<std::chrono::nanoseconds(Clock::now() - t1).count()/1e9<<endl;
//    }else{
//
//
//
//    }
//}

void RBF_Core::Set_HermiteRBF(vector<double>&pts){
#ifdef DEBUG
    cout<<"Set_HermiteRBF"<<endl;
#endif
    //for(auto a:pts)cout<<a<<' ';cout<<endl;
    isHermite = true;

    const int dim = point_dimension;

    a.set_size(npt*(dim+1));
    M.set_size(npt*(dim+1),npt*(dim+1));
    double *p_pts = pts.data();

    //M00
    for(int i=0;i<npt;++i){
        for(int j=i;j<npt;++j){
            M(i,j) = M(j,i) = Kernal_Function_2p(p_pts+i*dim, p_pts+j*dim);
        }
    }

    //if(User_Lamnbda!=0)for(int i=0;i<npt;++i)M(i,i) += User_Lamnbda;

    //M01
//    double G[3];
    //double G[dim];
    std::vector<double> G(dim);
    for(int i=0;i<npt;++i){
        for(int j=0;j<npt;++j){

            Kernal_Gradient_Function_2p(p_pts+i*dim, p_pts+j*dim, G.data());


            for(int k=0;k<dim;++k)M(i,npt+j+k*npt) = G[k];
            for(int k=0;k<dim;++k)M(npt+j+k*npt,i) = G[k];

        }
    }

    //M11
//    double H[9];
    //double H[dim*dim];
    std::vector<double> H(dim * dim);
    for(int i=0;i<npt;++i){
        for(int j=i;j<npt;++j){

            Kernal_Hessian_Function_2p(p_pts+i*dim, p_pts+j*dim, H.data());


            for(int k=0;k<dim;++k)
                for(int l=0;l<dim;++l)
                    M(npt+j+l*npt,npt+i+k*npt) = M(npt+i+k*npt,npt+j+l*npt) = -H[k*dim+l];
        }
    }

    //cout<<std::setprecision(5)<<std::fixed<<M<<endl;

    bsize= dim+1;
    N.zeros(npt*(dim+1),dim+1);
    b.set_size(dim+1);

    for(int i=0;i<npt;++i){
        N(i,0) = 1;
        for(int j=0;j<dim;++j)N(i,j+1) = pts[i*dim+j];
    }
    for(int i=0;i<npt;++i){
        for(int j=0;j<dim;++j)N(npt+i+j*npt,j+1) = -1;
    }

    //cout<<N<<endl;
    //arma::vec eigval = eig_sym( M ) ;
    //cout<<eigval.t()<<endl;


    if(!isnewformula){
        cout<<"start solve M: "<<endl;
        auto t1 = Clock::now();
        if(isinv)Minv = inv(M);
        else {
            arma::mat Eye;
            Eye.eye(npt*(dim+1),npt*(dim+1));
            Minv = solve(M,Eye);
        }
        cout<<"solved M: "<<(invM_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;

        t1 = Clock::now();
        if(isinv)bprey = inv_sympd(N.t() * Minv * N) * N.t() * Minv;
        else {
            arma::mat Eye2;
            Eye2.eye(bsize,bsize);
            bprey = solve(N.t() * Minv * N, Eye2) * N.t() * Minv;
        }
        cout<<"solved bprey "<<std::chrono::nanoseconds(Clock::now() - t1).count()/1e9<<endl;
    }else{

    }
}


double Gaussian_2p(const double *p1, const double *p2, double sigma){

    return exp(-MyUtility::vecSquareDist(p1,p2)/(2*sigma*sigma));
}



void RBF_Core::Set_Actual_User_LSCoef(double user_ls){

    User_Lamnbda = User_Lamnbda_inject = user_ls > 0 ?  user_ls : 0;

}

void RBF_Core::Set_Actual_Hermite_LSCoef(double hermite_ls){

    ls_coef = Hermite_ls_weight_inject = hermite_ls > 0?hermite_ls:0;
}

void RBF_Core::Set_SparsePara(double spa){
    sparse_para = spa;
}

void RBF_Core::Set_User_Lamnda_ToMatrix(double user_ls){


    {
        Set_Actual_User_LSCoef(user_ls);
        auto t1 = Clock::now();
#ifdef DEBUG
        cout<<"setting K, HermiteApprox_Lamnda"<<endl;
#endif
        if(User_Lamnbda>0){
            arma::sp_mat eye;
            eye.eye(npt,npt);

            dI = inv(eye + User_Lamnbda*K00);
            saveK_finalH = K = K11 - (User_Lamnbda)*(K01.t()*dI*K01);

        }else saveK_finalH = K = K11;
#ifdef DEBUG
        cout<<"solved: "<<(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;
#endif
    }

    finalH = saveK_finalH;

}

void RBF_Core::Set_HermiteApprox_Lamnda(double hermite_ls){


    {
        Set_Actual_Hermite_LSCoef(hermite_ls);
        auto t1 = Clock::now();
#ifdef DEBUG
        cout<<"setting K, HermiteApprox_Lamnda"<<endl;
#endif
        if(ls_coef>0){
            arma::sp_mat eye;
            eye.eye(npt,npt);

            if(ls_coef > 0){
                arma:: mat tmpdI = inv(eye + (ls_coef+User_Lamnbda)*K00);
                K = K11 - (ls_coef+User_Lamnbda)*(K01.t()*tmpdI*K01);
            }else{
                K = saveK_finalH;
            }
        }
#ifdef DEBUG
        cout<<"solved: "<<(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;
#endif
    }


}



void RBF_Core::Set_Hermite_PredictNormal(vector<double>&pts){
    Set_HermiteRBF(pts);

    auto t1 = Clock::now();
#ifdef DEBUG
    cout<<"setting K"<<endl;
#endif

    int dim = point_dimension;

    if(!isnewformula){
        arma::mat D = N.t()*Minv;
        K = Minv - D.t()*inv(D*N)*D;
//        K = K.submat( npt, npt, npt*4-1, npt*4-1 );
        K = K.submat( npt, npt, npt*(dim+1)-1, npt*(dim+1)-1 );
        finalH = saveK_finalH = K;

    }else{  // default branch
#ifdef DEBUG
        cout<<"using new formula"<<endl;
#endif
        // bigM: matrix A in the paper
//        bigM.zeros((npt+1)*4,(npt+1)*4);
//        bigM.submat(0,0,npt*4-1,npt*4-1) = M;
//        bigM.submat(0,npt*4,(npt)*4-1, (npt+1)*4-1) = N;
//        bigM.submat(npt*4,0,(npt+1)*4-1, (npt)*4-1) = N.t();
        bigM.zeros((npt+1)*(dim+1),(npt+1)*(dim+1));
        bigM.submat(0,0,npt*(dim+1)-1,npt*(dim+1)-1) = M;
        bigM.submat(0,npt*(dim+1),(npt)*(dim+1)-1, (npt+1)*(dim+1)-1) = N;
        bigM.submat(npt*(dim+1),0,(npt+1)*(dim+1)-1, (npt)*(dim+1)-1) = N.t();

        //for(int i=0;i<4;++i)bigM(i+(npt)*4,i+(npt)*4) = 1;

        auto t2 = Clock::now();
        bigMinv = inv(bigM);
        setK_time= std::chrono::nanoseconds(Clock::now() - t2).count()/1e9;
#ifdef DEBUG
        cout<<"bigMinv: "<<setK_time<<endl;
#endif
		bigM.clear();
		// Minv: matrix J in the paper
		// Ninv: matrix K in the paper
//        Minv = bigMinv.submat(0,0,npt*4-1,npt*4-1);
//        Ninv = bigMinv.submat(0,npt*4,(npt)*4-1, (npt+1)*4-1);
        Minv = bigMinv.submat(0,0,npt*(dim+1)-1,npt*(dim+1)-1);
        Ninv = bigMinv.submat(0,npt*(dim+1),(npt)*(dim+1)-1, (npt+1)*(dim+1)-1);

        bigMinv.clear();
        //K = Minv - Ninv *(N.t()*Minv);
        K = Minv;
        K00 = K.submat(0,0,npt-1,npt-1);
//        K01 = K.submat(0,npt,npt-1,npt*4-1);
//        K11 = K.submat( npt, npt, npt*4-1, npt*4-1 );
        K01 = K.submat(0,npt,npt-1,npt*(dim+1)-1);
        K11 = K.submat( npt, npt, npt*(dim+1)-1, npt*(dim+1)-1 );

        M.clear();N.clear();
#ifdef DEBUG
        cout<<"K11: "<<K11.n_cols<<endl;
#endif


        //Set_Hermite_DesignedCurve();

        Set_User_Lamnda_ToMatrix(User_Lamnbda_inject);

		
//		arma::vec eigval, ny;
//		arma::mat eigvec;
//		ny = eig_sym( eigval, eigvec, K);
//		cout<<ny<<endl;
#ifdef DEBUG
        cout<<"K: "<<K.n_cols<<endl;
#endif
    }

    //K = ( K.t() + K )/2;
#ifdef DEBUG
    cout<<"solve K total: "<<(setK_time= std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;
#endif
    return;

}



void RBF_Core::SetInitnormal_Uninorm(){

    initnormals_uninorm = initnormals;

    int dim = point_dimension;
    for(int i=0;i<npt;++i)MyUtility::normalize(initnormals_uninorm.data()+i*dim,dim);

}

int RBF_Core::Solve_Hermite_PredictNormal_UnitNorm(){

    arma::vec eigval, ny;
    arma::mat eigvec;

    if(!isuse_sparse){
        ny = eig_sym( eigval, eigvec, K);
    }else{
//		cout<<"use sparse eigen"<<endl;
//        int k = 4;
//        do{
//            ny = eigs_sym( eigval, eigvec, sp_K, k, "sa" );
//            k+=4;
//        }while(ny(0)==0);
    }

#ifdef DEBUG
    cout<<"eigval(0): "<<eigval(0)<<endl;
#endif

    int smalleig = 0;

    int dim = point_dimension;

    initnormals.resize(npt*dim);
    arma::vec y(npt*(dim+1));
    for(int i=0;i<npt;++i)y(i) = 0;
    for(int i=0;i<npt*dim;++i)y(i+npt) = eigvec(i,smalleig);
    for(int i=0;i<npt;++i){
        for (int j=0; j<dim; ++j) {
            initnormals[i*dim+j] = y(npt+i+npt*j);
        }
    }

    SetInitnormal_Uninorm();
#ifdef DEBUG
    cout<<"Solve_Hermite_PredictNormal_UnitNorm finish"<<endl;
#endif
    return 1;
}



/***************************************************************************************************/
/***************************************************************************************************/
double acc_time;

static int countopt = 0;
double optfunc_Hermite(const vector<double>&x, vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    RBF_Core *drbf = reinterpret_cast<RBF_Core*>(fdata);
    int n = drbf->npt;
    arma::vec arma_x(n*3);

    //(  sin(a)cos(b), sin(a)sin(b), cos(a)  )  a =>[0, pi], b => [-pi, pi];
    vector<double>sina_cosa_sinb_cosb(n * 4);
    for(int i=0;i<n;++i){
        int ind = i*4;
        sina_cosa_sinb_cosb[ind] = sin(x[i*2]);
        sina_cosa_sinb_cosb[ind+1] = cos(x[i*2]);
        sina_cosa_sinb_cosb[ind+2] = sin(x[i*2+1]);
        sina_cosa_sinb_cosb[ind+3] = cos(x[i*2+1]);
    }

    for(int i=0;i<n;++i){
        auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        //        int ind = i*3;
        //        arma_x(ind) = p_scsc[0] * p_scsc[3];
        //        arma_x(ind+1) = p_scsc[0] * p_scsc[2];
        //        arma_x(ind+2) = p_scsc[1];
        arma_x(i) = p_scsc[0] * p_scsc[3];
        arma_x(i+n) = p_scsc[0] * p_scsc[2];
        arma_x(i+n*2) = p_scsc[1];
    }

    arma::vec a2;
    //if(drbf->isuse_sparse)a2 = drbf->sp_H * arma_x;
    //else
    a2 = drbf->finalH * arma_x;


    if (!grad.empty()) {

        grad.resize(n*2);

        for(int i=0;i<n;++i){
            auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;

            //            int ind = i*3;
            //            grad[i*2] = a2(ind) * p_scsc[1] * p_scsc[3] + a2(ind+1) * p_scsc[1] * p_scsc[2] - a2(ind+2) * p_scsc[0];
            //            grad[i*2+1] = -a2(ind) * p_scsc[0] * p_scsc[2] + a2(ind+1) * p_scsc[0] * p_scsc[3];

            grad[i*2] = a2(i) * p_scsc[1] * p_scsc[3] + a2(i+n) * p_scsc[1] * p_scsc[2] - a2(i+n*2) * p_scsc[0];
            grad[i*2+1] = -a2(i) * p_scsc[0] * p_scsc[2] + a2(i+n) * p_scsc[0] * p_scsc[3];

        }
    }

    double re = arma::dot( arma_x, a2 );
    countopt++;

    acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);

    //cout<<countopt++<<' '<<re<<endl;
    return re;

}

double optfunc_Hermite_2D(const vector<double>&x, vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    RBF_Core *drbf = reinterpret_cast<RBF_Core*>(fdata);
    int n = drbf->npt;
    int dim = 2;
    arma::vec arma_x(n*dim);

    //(  cos(b), sin(b)  )  b => [-pi, pi];
    vector<double> sinb_cosb(n * 2);
    for(int i=0;i<n;++i){
        int ind = i*2;
        sinb_cosb[ind] = sin(x[i]);
        sinb_cosb[ind+1] = cos(x[i]);
    }

    for(int i=0;i<n;++i){
        auto p_scsc = sinb_cosb.data()+i*2;
        //
        arma_x(i) = p_scsc[1];
        arma_x(i+n) = p_scsc[0];
    }

    arma::vec a2;
    //if(drbf->isuse_sparse)a2 = drbf->sp_H * arma_x;
    //else
    a2 = drbf->finalH * arma_x;


    if (!grad.empty()) {

        grad.resize(n);

        for(int i=0;i<n;++i){
            auto p_scsc = sinb_cosb.data()+i*2;

            grad[i] = -a2(i) * p_scsc[0] + a2(i+n) * p_scsc[1];
        }
    }

    double re = arma::dot( arma_x, a2 );
    countopt++;

    acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);

    //cout<<countopt++<<' '<<re<<endl;
    return re;

}



int RBF_Core::Opt_Hermite_PredictNormal_UnitNormal(){


    sol.solveval.resize(npt * 2);

    for(int i=0;i<npt;++i){
        double *veccc = initnormals.data()+i*3;
        {
            //MyUtility::normalize(veccc);
            sol.solveval[i*2] = atan2(sqrt(veccc[0]*veccc[0]+veccc[1]*veccc[1]),veccc[2] );
            sol.solveval[i*2 + 1] = atan2( veccc[1], veccc[0]   );
        }

    }
    //cout<<"smallvec: "<<smallvec<<endl;

    if(1){
        vector<double>upper(npt*2);
        vector<double>lower(npt*2);
        for(int i=0;i<npt;++i){
            upper[i*2] = 2 * my_PI;
            upper[i*2 + 1] = 2 * my_PI;

            lower[i*2] = -2 * my_PI;
            lower[i*2 + 1] = -2 * my_PI;
        }

        countopt = 0;
        acc_time = 0;

        //LocalIterativeSolver(sol,kk==0?normals:newnormals,300,1e-7);
        Solver::nloptwrapper(lower,upper,optfunc_Hermite,this,1e-7,3000,sol);
#ifdef DEBUG
        cout<<"number of call: "<<countopt<<" t: "<<acc_time<<" ave: "<<acc_time/countopt<<endl;
#endif
        callfunc_time = acc_time;
        solve_time = sol.time;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;

    }
    newnormals.resize(npt*3);
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;
    for(int i=0;i<npt;++i){

        double a = sol.solveval[i*2], b = sol.solveval[i*2+1];
        newnormals[i*3]   = y(npt+i) = sin(a) * cos(b);
        newnormals[i*3+1] = y(npt+i+npt) = sin(a) * sin(b);
        newnormals[i*3+2] = y(npt+i+npt*2) = cos(a);
        MyUtility::normalize(newnormals.data()+i*3);
    }

    Set_RBFCoef(y);

    //sol.energy = arma::dot(a,M*a);
#ifdef DEBUG
    cout<<"Opt_Hermite_PredictNormal_UnitNormal"<<endl;
#endif
    return 1;
}

int RBF_Core::Opt_Hermite_PredictNormal_UnitNormal_2D(){

    sol.solveval.resize(npt);

    int dim = point_dimension;

    for(int i=0;i<npt;++i){
        double *veccc = initnormals.data()+i*dim;
        {
            sol.solveval[i] = atan2( veccc[1], veccc[0]);
        }
    }

    if(1){
        vector<double>upper(npt);
        vector<double>lower(npt);
        for(int i=0;i<npt;++i){
            upper[i] = 2 * my_PI;
            lower[i] = -2 * my_PI;
        }

        countopt = 0;
        acc_time = 0;

        //LocalIterativeSolver(sol,kk==0?normals:newnormals,300,1e-7);
        Solver::nloptwrapper(lower,upper,optfunc_Hermite_2D,this,1e-7,3000,sol);
        cout<<"number of call: "<<countopt<<" t: "<<acc_time<<" ave: "<<acc_time/countopt<<endl;
        callfunc_time = acc_time;
        solve_time = sol.time;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;

    }
    newnormals.resize(npt*dim);
    arma::vec y(npt*(dim+1));
    for(int i=0;i<npt;++i)y(i) = 0;
    for(int i=0;i<npt;++i){
        double b = sol.solveval[i];
        newnormals[i*dim] = y(npt+i) = cos(b);
        newnormals[i*dim+1] = y(npt+i+npt) = sin(b);
        MyUtility::normalize(newnormals.data()+i*dim, dim);
    }

    Set_RBFCoef(y);

    //sol.energy = arma::dot(a,M*a);
    cout<<"Opt_Hermite_PredictNormal_UnitNormal"<<endl;
    return 1;
}

void RBF_Core::Set_RBFCoef(arma::vec &y){
#ifdef DEBUG
    cout<<"Set_RBFCoef"<<endl;
#endif
    if(curMethod==HandCraft){
        cout<<"HandCraft, not RBF"<<endl;
        return;
    }
    if(!isnewformula){
        b = bprey * y;
        a = Minv * (y - N*b);
    }else{
        int dim = point_dimension;
        if(User_Lamnbda>0)y.subvec(0,npt-1) = -User_Lamnbda*dI*K01*y.subvec(npt,npt*(dim+1)-1);

        a = Minv*y;
        b = Ninv.t()*y;

    }
}

bool RBF_Core::Write_RBFCoeff(string fname)
{
    ofstream outer(fname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output file " << fname << endl;
        return false;
    }

    for (int i = 0; i < a.size(); ++i) {
        outer << a(i) << " ";
    }
    outer << endl;

    for (int i = 0; i < b.size(); ++i) {
        outer << b(i) << " ";
    }
    outer << endl;

    outer.close();
    cout<<"saving finish: "<< fname <<endl;
    return true;

}

bool RBF_Core::Write_tmp(string fname)
{
    ofstream outer(fname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output file " << fname << endl;
        return false;
    }

    // Minv
    for (int i = 0; i < Minv.n_rows; ++i) {
        for (int j=0; j < Minv.n_cols; ++j) {
            outer << Minv(i,j) << " ";
        }
        outer << endl;
    }
    outer << endl;

    // Ninv
    for (int i = 0; i < Ninv.n_rows; ++i) {
        for (int j=0; j < Ninv.n_cols; ++j) {
            outer << Ninv(i,j) << " ";
        }
        outer << endl;
    }
    outer << endl;

    outer.close();
    cout<<"saving finish: "<< fname <<endl;
    return true;

}


int RBF_Core::Lamnbda_Search_GlobalEigen(){

    vector<double>lamnbda_list({0, 0.001, 0.01, 0.1, 1});
    vector<double>initen_list(lamnbda_list.size());
    vector<double>finalen_list(lamnbda_list.size());
    vector<vector<double>>init_normallist;
    vector<vector<double>>opt_normallist;

    lamnbda_list_sa = lamnbda_list;
    for(int i=0;i<lamnbda_list.size();++i){

        Set_HermiteApprox_Lamnda(lamnbda_list[i]);

        if(curMethod==Hermite_UnitNormal){
            Solve_Hermite_PredictNormal_UnitNorm();
        }

        //Solve_Hermite_PredictNormal_UnitNorm();
        OptNormal(1);

        initen_list[i] = sol.init_energy;
        finalen_list[i] = sol.energy;

        init_normallist.emplace_back(initnormals);
        opt_normallist.emplace_back(newnormals);
    }

    lamnbdaGlobal_Be.emplace_back(initen_list);
    lamnbdaGlobal_Ed.emplace_back(finalen_list);

    cout<<std::setprecision(8);
#ifdef DEBUG
    for(int i=0;i<initen_list.size();++i){
        cout<<lamnbda_list[i]<<": "<<initen_list[i]<<" -> "<<finalen_list[i]<<endl;
    }
#endif

    int minind = min_element(finalen_list.begin(),finalen_list.end()) - finalen_list.begin();
#ifdef DEBUG
    cout<<"min energy: "<<endl;
    cout<<lamnbda_list[minind]<<": "<<initen_list[minind]<<" -> "<<finalen_list[minind]<<endl;
#endif


    initnormals = init_normallist[minind];
    SetInitnormal_Uninorm();
    newnormals = opt_normallist[minind];
	return 1;
}




void RBF_Core::Print_LamnbdaSearchTest(string fname){


    cout<<setprecision(7);
    cout<<"Print_LamnbdaSearchTest"<<endl;
    for(int i=0;i<lamnbda_list_sa.size();++i)cout<<lamnbda_list_sa[i]<<' ';cout<<endl;
    cout<<lamnbdaGlobal_Be.size()<<endl;
    for(int i=0;i<lamnbdaGlobal_Be.size();++i){
        for(int j=0;j<lamnbdaGlobal_Be[i].size();++j){
            cout<<lamnbdaGlobal_Be[i][j]<<"\t"<<lamnbdaGlobal_Ed[i][j]<<"\t";
        }
        cout<<gtBe[i]<<"\t"<<gtEd[i]<<endl;
    }

    ofstream fout(fname);
    fout<<setprecision(7);
    if(!fout.fail()){
        for(int i=0;i<lamnbda_list_sa.size();++i)fout<<lamnbda_list_sa[i]<<' ';fout<<endl;
        fout<<lamnbdaGlobal_Be.size()<<endl;
        for(int i=0;i<lamnbdaGlobal_Be.size();++i){
            for(int j=0;j<lamnbdaGlobal_Be[i].size();++j){
                fout<<lamnbdaGlobal_Be[i][j]<<"\t"<<lamnbdaGlobal_Ed[i][j]<<"\t";
            }
            fout<<gtBe[i]<<"\t"<<gtEd[i]<<endl;
        }
    }
    fout.close();

}


