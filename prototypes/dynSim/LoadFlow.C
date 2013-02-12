


class LoadFlow {
	Mat ymat;
	Vec vthera;
	Vec pq;

	void solve();
	void yswitch();
	void 
};


LoadFlow::LoadFlow(Mat& ymat) {
	
	
};

void LoadFlow::solve(Vec& pq,Vec& vtheta)
{

	KSPConvergedReason reason;
	KSP ksp;

	KSPCreate(&ksp);
	KSPSetOperators(ksp,Y,Y,DIFFERENT_NONZERO_PATTERN);
//      KSPSetType(ksp,KSPCG);
	KSPSetTolerance(ksp,eps,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);
	KSPSetFromOptions(ksp);
	KSPSolve(ksp,pq,vtheta);
	KSPGetConvergedReason(ksp,&reason);
	if (reason < 0) {
		cerr << " did not converge\n";
		utils::gridpack_exit();
	}else{
		iter=KSPGetIterations(ksp,&its);
		cerr << "converged in " << iter << " iterations\n";
	}	
	KSPDestroy(ksp);
}