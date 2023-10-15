#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>

using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS

void invTransformH( vec eigval, mat &H ){
	const size_t num = sum( eigval < 1e-8 );
	const uvec idx = find( eigval < 1e-8 );
	for(size_t i=0; i<num; i++){ eigval[idx[i]] = 1e-5; }
	eigval += 0.5*randu( H.n_rows );
	mat Q(H.n_rows,H.n_cols), R(H.n_rows,H.n_cols), O(H.n_rows, H.n_cols);
	qr( Q, R, H );
	O = Q * diagmat( R.diag()/abs(R.diag()) );
	H = O.t() * diagmat(1.0/eigval) * O;
	return;
}

//*************************************************************************************//
//                             AVERAGE INFORMATION METHOD                              //
//*************************************************************************************//

//***// [[Rcpp::export]]是必须的
// [[Rcpp::export]]
SEXP glmmAI(SEXP Yin, SEXP Xin, SEXP numKin, SEXP Phiin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin)
{/*Average Information*/
	try {
		vec Y = as<vec>(Yin);
		mat X = as<mat>(Xin);
		int numK = Rcpp::as<int>(numKin);
		const Rcpp::List Phi(Phiin);
		vec D = as<vec>(Din);//这里的D实际上是D^2
		vec tau = as<vec>(tauin);
		const uvec fixtau = as<uvec>(fixtauin);
		int numK2 = sum( fixtau == 0 );
		const double tol = Rcpp::as<double>(tolin);
		uvec ZERO = (tau < tol);
		size_t numIDV = X.n_rows, numCVT = X.n_cols;
		mat Hinv(numCVT, numCVT), XtHinvX_inv(numCVT, numCVT), P(numIDV, numIDV);
		vec alpha(numCVT), eta(numIDV);

		mat H = tau[0] * diagmat(1.0/D);// Phi里面的第一个为 I 或 W 或 D 矩阵
		for (size_t i=1; i<=numK; ++i) {
			H = H + tau[i] * as<mat>(Phi[i-1]);//第一个表示残差；第二个表示遗传方差
		}

		Hinv = inv_sympd (H);
		mat HinvX = Hinv * X;
		mat XtHinvX = X.t() * HinvX;
		XtHinvX_inv = inv_sympd (XtHinvX);
		P = Hinv - HinvX * XtHinvX_inv * HinvX.t();
		alpha = XtHinvX_inv * HinvX.t() * Y;
		eta = Y - tau[0] * (Hinv * (Y - X * alpha)) / D;

		mat K0 = diagmat(1.0 / D);
		mat FisherMat;
		if (numK2 > 0) {
			const uvec idxtau = find(fixtau == 0);
			mat AImat(numK2, numK2);//average information matrix
			vec PY = P * Y;
			vec score(numK2), PKPY;
			//mat PK0 = P * K0;//NEW
			double tr1 = sum(diagvec(P) / D);
			for (size_t i=0; i<numK2; ++i) {
				if (i == 0 && idxtau[0] == 0) {
					PKPY = PY / D;
					score[0]    = 0.5 * (dot(PKPY, PY) - tr1);
					AImat(0, 0) = 0.5 * dot(PKPY, P * PKPY); //AI
					//AImat(0, 0) = -0.5 * accu(PK0 % PK0) + dot(PKPY, P * PKPY);//OI
					//AImat(0, 0) =  0.5 * accu(PK0 % PK0);//EI
				} else {
				PKPY = P * as<mat>(Phi[idxtau[i]-1]) * PY;
					score[i] = 0.5 * (dot(Y, PKPY) - accu(P % as<mat>(Phi[idxtau[i]-1])));
					for(size_t j=0; j<=i; ++j) {
						if(j == 0 && idxtau[0] == 0) {
							AImat(i, 0) = 0.5 * dot(PY / D, PKPY);//AI
							//AImat(i, 0) = -0.5 * accu(PK0 % (P * as<mat>(Phi[idxtau[i]-1]))) + dot(PY / D, PKPY);//OI
							//AImat(i, 0) =  0.5 * accu(PK0 % (P * as<mat>(Phi[idxtau[i]-1])));//EI
							AImat(0, i) = AImat(i, 0);
						} else {
							AImat(i, j) = 0.5 * dot(PY, as<mat>(Phi[idxtau[j]-1]) * PKPY);//AI
							//AImat(i, j) = -0.5 * accu((P * as<mat>(Phi[idxtau[j]-1])) % (P * as<mat>(Phi[idxtau[i]-1]))) + dot(PY, as<mat>(Phi[idxtau[j]-1]) * PKPY);//OI
							//AImat(i, j) = 0.5 * accu((P * as<mat>(Phi[idxtau[j]-1])) % (P * as<mat>(Phi[idxtau[i]-1])));//EI
							if(j!=i) {AImat(j, i) = AImat(i, j);}
						}
					}//end for j
				}
			}// end for i
			vec Dtau = solve(AImat, score);
			vec tau0 = tau;
			tau.elem( idxtau ) = tau0.elem( idxtau ) + Dtau;
			tau.elem( find(ZERO % (tau < tol)) ).zeros();
			double step = 1.0;
			while(any(tau < 0.0)) {
			    step *= 0.5;
				tau.elem( idxtau ) = tau0.elem( idxtau ) + step * Dtau;
				tau.elem( find(ZERO % (tau < tol)) ).zeros();
			}
			tau.elem( find(tau < tol) ).zeros();
			FisherMat = AImat;
		}
		// return values
		return List::create(Named("tau") = tau, Named("P") = P, Named("cov") = XtHinvX_inv, 
			Named("alpha") = alpha, Named("eta") = eta, 
			Named("H")=H, Named("FisherMat")=FisherMat);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}
/////////////////////////////////////////////////////////////////////////////////////////
//                             CODE END HERE                                           //
/////////////////////////////////////////////////////////////////////////////////////////