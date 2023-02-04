#include "PIC.h"

void PIC_TYP::MPI_AllreduceVec(const params_TYP * params, arma::vec * v)
{
	arma::vec recvbuf = zeros(v->n_elem);

	MPI_Allreduce(v->memptr(), recvbuf.memptr(), v->n_elem, MPI_DOUBLE, MPI_SUM, params->mpi.MPI_TOPO);

	*v = recvbuf;
}

void PIC_TYP::MPI_SendVec(const params_TYP * params, arma::vec * v)
{
	// Send vector from PARTICLE ROOT to FIELDS ROOT:
	if (params->mpi.IS_PARTICLES_ROOT)
    {
		MPI_Send(v->memptr(), v->n_elem, MPI_DOUBLE, params->mpi.FIELDS_ROOT_WORLD_RANK, PARTICLES_TAG, MPI_COMM_WORLD);
	}

	if (params->mpi.IS_FIELDS_ROOT)
    {
		MPI_Recv(v->memptr(), v->n_elem, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, PARTICLES_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// Broadcast from FIELDS ROOT to all other FIELDS MPIs:
	if (params->mpi.COMM_COLOR == FIELDS_MPI_COLOR)
    {
        MPI_Bcast(v->memptr(), v->n_elem, MPI_DOUBLE, 0, params->mpi.COMM);
    }
}


void PIC_TYP::MPI_ReduceVec(const params_TYP * params, arma::vec * v)
{
    // Create receive buffer:
    // ======================
    arma::vec recvbuf = zeros(v->n_elem);

    // Reduce the vector at the PARTICLE ROOT:
    // ======================================
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        MPI_Reduce(v->memptr(), recvbuf.memptr(), v->n_elem, MPI_DOUBLE, MPI_SUM, 0, params->mpi.COMM);

        if (params->mpi.IS_PARTICLES_ROOT)
        {
             *v = recvbuf;
        }
    }
}

void PIC_TYP::MPI_Allgathervec(const params_TYP * params, arma::vec * field)
{
	unsigned int iIndex = params->mpi.iIndex;
	unsigned int fIndex = params->mpi.fIndex;

	arma::vec recvbuf(params->mesh_params.Nx);
	arma::vec sendbuf(params->mesh_params.Nx_PER_MPI);

	sendbuf = field->subvec(iIndex, fIndex);
	MPI_Allgather(sendbuf.memptr(), params->mesh_params.Nx_PER_MPI, MPI_DOUBLE, recvbuf.memptr(), params->mesh_params.Nx_PER_MPI, MPI_DOUBLE, params->mpi.MPI_TOPO);
	field->subvec(1, params->mesh_params.Nx) = recvbuf;
}


void PIC_TYP::MPI_Recvvec(const params_TYP * params, arma::vec * field)
{
	// We send the vector from root process of fields to root process of particles
	arma::vec recvbuf(params->mesh_params.Nx);
	arma::vec sendbuf(params->mesh_params.Nx);

	sendbuf = field->subvec(1, params->mesh_params.Nx);

	if (params->mpi.IS_FIELDS_ROOT)
    {
		MPI_Send(sendbuf.memptr(), params->mesh_params.Nx, MPI_DOUBLE, params->mpi.PARTICLES_ROOT_WORLD_RANK, 0, MPI_COMM_WORLD);
	}

	if (params->mpi.IS_PARTICLES_ROOT)
    {
		MPI_Recv(recvbuf.memptr(), params->mesh_params.Nx, MPI_DOUBLE, params->mpi.FIELDS_ROOT_WORLD_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		field->subvec(1, params->mesh_params.Nx) = recvbuf;
	}

	// Then, the fields is broadcasted to all processes in the particles communicator COMM
	if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
		sendbuf = field->subvec(1, params->mesh_params.Nx);

		MPI_Bcast(sendbuf.memptr(), params->mesh_params.Nx, MPI_DOUBLE, 0, params->mpi.COMM);

		field->subvec(1, params->mesh_params.Nx) = sendbuf;
	}
}

void PIC_TYP::MPI_Recv_AllFields(const params_TYP * params, fields_TYP * fields)
{
	// Send field data from FIELDS ranks and recieve fields data at PARTICLE ranks
	MPI_Recvvec(params,&fields->Ex_m);
	MPI_Recvvec(params,&fields->Bx_m);
	MPI_Recvvec(params,&fields->dBx_m);

	if (params->SW.RFheating == 1)
	{
		MPI_Recvvec(params,&fields->ddBx_m);
	}
}

void PIC_TYP::fillGhosts(arma::vec * C)
{
	int Nx = C->n_elem;

	(*C)(0)    = (*C)(1);
	(*C)(Nx-1) = (*C)(Nx-2);
}

void PIC_TYP::fillGhost_AllFields(const params_TYP * params, fields_TYP * fields)
{
	fillGhosts(&fields->Ex_m);
	fillGhosts(&fields->Bx_m);
	fillGhosts(&fields->dBx_m);

	if (params->SW.RFheating == 1)
	{
		fillGhosts(&fields->ddBx_m);
	}
}

void PIC_TYP::fill4Ghosts(arma::vec * v)
{
	int Nx = v->n_elem;

	v->subvec(0,1)       = v->subvec(2,3);
	v->subvec(Nx-2,Nx-1) = v->subvec(Nx-4,Nx-3);
}

void PIC_TYP::smooth(arma::vec * v, double as)
{
	int Nx(v->n_elem);

	arma::vec b = zeros(Nx);

	double wc(0.75); 	// center weight
	double ws(0.125);	// sides weight

	//Step 1: Averaging process
	b.subvec(1, Nx-2) = v->subvec(1, Nx-2);

	fillGhosts(&b);

	b.subvec(1, Nx-2) = wc*b.subvec(1, Nx-2) + ws*b.subvec(2, Nx-1) + ws*b.subvec(0, Nx-3);

	//Step 2: Averaged weighted variable estimation.
	v->subvec(1, Nx-2) = (1.0 - as)*v->subvec(1, Nx-2) + as*b.subvec(1, Nx-2);
}

// Constructor:
PIC_TYP::PIC_TYP(const params_TYP * params, const mesh_TYP * mesh, fields_TYP * fields, vector<ions_TYP> * IONS, electrons_TYP * electrons)
{
	// Get latest mesh-defined values from FIELDS ranks:
	// =================================================
	MPI_Recv_AllFields(params,fields);

	// Fill the ghost cells in all fields:
	fillGhost_AllFields(params,fields);

	// Assign cell for all particles:
  // =============================
	assignCell_AllSpecies(params,mesh,IONS);

	// Interpolate all fields on all species:
	// ======================================
	interpolateFields_AllSpecies(params,IONS,fields);

	// Interpolate electron temperature on all species:
	// ===============================================
	interpolateElectrons_AllSpecies(params,IONS,electrons);

	// Calculate ion moments and populate mesh-defined ion moments:
	// ============================================================
	// Run 3 dummy cycles to load "n" and "nv" at previous time steps:
	for(int tt=0; tt<3; tt++)
	{
		extrapolateMoments_AllSpecies(params,fields,IONS);
	}
}

void PIC_TYP::interpolateScalarField(const params_TYP * params, ions_TYP * IONS, const arma::vec * F_m, arma::vec * F_p)
{
    int Nx =  params->mesh_params.Nx + 4; //Mesh size along the x axis (considering the gosht cell)
    int N_CP(IONS->N_CP_MPI);

    // Allocate memory and initialize to zero:
    arma::vec F = zeros(Nx);

    // Fill in vector F with mesh-defined data:
    F.subvec(1,Nx-2) = *F_m;

    // Take care of ghost cells:
    fill4Ghosts(&F);

    #pragma omp parallel for default(none) shared(params, IONS, F_p, F) firstprivate(N_CP)
    for(int ii=0; ii<N_CP; ii++)
    {
        int ix = IONS->mn(ii) + 2;

        (*F_p)(ii) += IONS->wxl(ii)*F(ix-1);
        (*F_p)(ii) += IONS->wxc(ii)*F(ix);
        (*F_p)(ii) += IONS->wxr(ii)*F(ix+1);

    }// omp parallel for
}

void PIC_TYP::interpolateFields(const params_TYP * params, ions_TYP * IONS, const fields_TYP * fields)
{
	// Need to use SW in order to enable/disable ddBX interpolation

	// Allocate memory for field variables:
	arma::vec Ex_p   = zeros(IONS->N_CP_MPI, 1);
	arma::vec Bx_p   = zeros(IONS->N_CP_MPI, 1);
	arma::vec dBx_p  = zeros(IONS->N_CP_MPI, 1);

	// Interpolate mesh-defined fields into particles:
	interpolateScalarField(params, IONS, &fields->Ex_m , &Ex_p );
	interpolateScalarField(params, IONS, &fields->Bx_m , &Bx_p );
	interpolateScalarField(params, IONS, &fields->dBx_m, &dBx_p);

	// Assign values:
	IONS->Ex_p  = Ex_p;
	IONS->Bx_p  = Bx_p;
	IONS->dBx_p = dBx_p;

	if (params->SW.RFheating == 1)
	{
		arma::vec ddBx_p = zeros(IONS->N_CP_MPI, 1);
		interpolateScalarField(params, IONS, &fields->ddBx_m, &ddBx_p);
		IONS->ddBx_p = ddBx_p;
	}
}

void PIC_TYP::interpolateFields_AllSpecies(const params_TYP * params, vector<ions_TYP> * IONS, const fields_TYP * fields)
{
	// Interpolate all fields on all species:
	// =======================
	for(int ss=0;ss<IONS->size();ss++)
	{
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
			// Interpolate mesh-defined fields into ALL particles locations:
			interpolateFields(params, &IONS->at(ss), fields);
		}
	}
}

void PIC_TYP::interpolateElectrons(const params_TYP * params, ions_TYP * IONS, const electrons_TYP * electrons)
{
	// Allocate memory:
	arma::vec Te_p = zeros(IONS->N_CP_MPI, 1);

	// Interpolate mesh-defined fields into particles:
	interpolateScalarField(params, IONS, &electrons->Te_m , &Te_p );

	// Assign values:
	IONS->Te_p = Te_p;
}

void PIC_TYP::interpolateElectrons_AllSpecies(const params_TYP * params, vector<ions_TYP> * IONS, const electrons_TYP * electrons)
{
	// Interpolate all fields on all species:
	// =======================
	for(int ss=0;ss<IONS->size();ss++)
	{
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
			// Interpolate mesh-defined electron temperature into ALL particles locations:
			interpolateElectrons(params, &IONS->at(ss), electrons);
		}
	}
}

void PIC_TYP::interpEM(const params_TYP * params, mesh_TYP * mesh, const ions_TYP * IONS, const fields_TYP * fields, arma::rowvec * ZN, arma::rowvec * EM)
{
    // Assign cell:
    double xp   = (*ZN)(0);
    double xMin = params->mesh_params.Lx_min;
    double dx   = params->mesh_params.dx;
    int m       = round( 0.5 + (xp - xMin)/dx ) - 1;

    // During RK4, some projections can be go out of bound
    if ( m >= params->mesh_params.Nx)
    {
        m = params->mesh_params.Nx - 1;
    }
    if ( m < 0)
    {
        m = 0;
    }

    // Distance to nearest grid point:
    double x = mesh->xm(m) - xp;

    // Assignment function:
    arma::vec W = zeros<vec>(3);
    W(0) = 0.5*pow( 1.5 + ((x - dx)/dx) ,2); // Left:
    W(1) = 0.75 - pow(x/dx,2);               // Center:
    W(2) = 0.5*pow(1.5 - ((x + dx)/dx) ,2 ); // Right:

    // Nearest grid point:
    int ix = m + 1;

    // Temporary storage for interpolated fields:
    arma::vec f = zeros<vec>(3);

    // Interpolate:
    // EX:
    f(0) = fields->Ex_m(ix - 1);
    f(1) = fields->Ex_m(ix);
    f(2) = fields->Ex_m(ix + 1);
    (*EM)(0) = arma::dot(f,W);

    // BX:
    f(0) = fields->Bx_m(ix - 1);
    f(1) = fields->Bx_m(ix);
    f(2) = fields->Bx_m(ix + 1);
    (*EM)(1) = arma::dot(f,W);

    // dBX:
    f(0) = fields->dBx_m(ix - 1);
    f(1) = fields->dBx_m(ix);
    f(2) = fields->dBx_m(ix + 1);
    (*EM)(2) = arma::dot(f,W);
}

void PIC_TYP::calculateF(const params_TYP * params, const ions_TYP * IONS, arma::rowvec * ZN, arma::rowvec * EM, arma::rowvec * F)
{
    // Ion parameters:
    double qa = IONS->Q;
    double Ma = IONS->M;

	// Gather fields:
    double E    = (*EM)(0);
    double B    = (*EM)(1);
    double dB   = (*EM)(2);

    // Output:
	if (params->advanceParticleMethod == 1)
	{
		// Gather particle states:
		double vpar = (*ZN)(1);
		double vper = (*ZN)(2);

		// Output:
	  (*F)(0) = +vpar;
		(*F)(1) = -0.5*vper*vper*dB/B + (qa/Ma)*E;
    (*F)(2) = +0.5*vper*vpar*dB/B;
	}
	if (params->advanceParticleMethod == 2)
	{
		// Gather particle states:
		double vpar = (*ZN)(1);
		double mu = (*ZN)(2);

	  (*F)(0) = +vpar;
		(*F)(1) = -(mu/Ma)*dB + (qa/Ma)*E;
    (*F)(2) = 0;
	}
}

void PIC_TYP::advanceParticles(const params_TYP * params, mesh_TYP * mesh, fields_TYP * fields, vector<ions_TYP> * IONS)
{
  // Get latest mesh-defined values from FIELDS MPIs:
	MPI_Recv_AllFields(params,fields);

	// Fill the ghost cells in all fields:
	fillGhost_AllFields(params,fields);

	// Iterate over all the ion species:
	for(int ss=0; ss<IONS->size(); ss++)
	{
		if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
		{
			// Number of particles:
			int N_CP = IONS->at(ss).N_CP_MPI;

      // Time step:
      double DT = params->DT;

			// Ion mass:
			double Ma = IONS->at(ss).M;

			#pragma omp parallel for default(none) shared(IONS, params, fields, DT, ss, Ma, std::cout, mesh) firstprivate(N_CP, F_C_DS)
      for(int ii=0;ii<N_CP;ii++)
			{
        // Start RK4 solution:
        //==============================================================

        // Equations of motion are expressed as follows:
        //
        // dZ/dt = F, thus Z1 = Z0 + dZ
				//
        // where
				//
        // dZ = F(Z0)*dt
        // Z = [z, vz, vp]
        // F(1) = +vz
        // F(2) = -0.5*vp*vp*dB/B + (q/Ma)*E
        // F(3) = +0.5*vp*vz*dB/B

				// Assemble vectors to use in RK4 method:
				arma::rowvec Z0  = zeros<rowvec>(3);
        arma::rowvec ZN  = zeros<rowvec>(3);
				arma::rowvec Z1  = zeros<rowvec>(3);
				arma::rowvec EM  = zeros<rowvec>(3);
        arma::rowvec F   = zeros<rowvec>(3);
        arma::rowvec dZ1 = zeros<rowvec>(3);
        arma::rowvec dZ2 = zeros<rowvec>(3);
        arma::rowvec dZ3 = zeros<rowvec>(3);
        arma::rowvec dZ4 = zeros<rowvec>(3);

				// Extract particle states:
				double x    = IONS->at(ss).x_p(ii);
				double vpar = IONS->at(ss).v_p(ii,0);
				double vper = IONS->at(ss).v_p(ii,1);

				// Extract particle-defined fields:
				double E  = IONS->at(ss).Ex_p(ii);
				double B  = IONS->at(ss).Bx_p(ii);
				double dB = IONS->at(ss).dBx_p(ii);

				// Initialize initial particle state Z0:
				Z0 = {x, vpar, vper};

				if ( isnan(Z0(0)) || isnan(Z0(1)) || isnan(Z0(2)) )
				{
					cout << "*Z0(0), 1 = " << Z0(0) << endl;
					cout << "*Z0(1), 1 = " << Z0(1) << endl;
					cout << "*Z0(2), 1 = " << Z0(2) << endl;
				}

				// Interpolate fields at current particle postion:
				interpEM(params, mesh, &IONS->at(ss), fields, &Z0, &EM);

				// Select solution method:
				switch (params->advanceParticleMethod)
				{
				case 1:
					// Do nothing
					break;
				case 2:
					// Use magnetic moment
					B  = EM(1);
					double mu = 0.5*Ma*pow(vper,2)/B;
					Z0(2) = mu;
					break;
				}

        // Step 1:
        ZN = Z0;
        calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
        dZ1 = F*DT;

        // Step 2:
        ZN = Z0 + dZ1/2;
        interpEM(params, mesh, &IONS->at(ss), fields, &ZN, &EM);
        calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
        dZ2 = F*DT;

        // Step 3:
        ZN = Z0 + dZ2/2;
        interpEM(params, mesh, &IONS->at(ss), fields, &ZN, &EM);
        calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
        dZ3 = F*DT;

        // Step 4:
        ZN = Z0 + dZ3;
        interpEM(params, mesh, &IONS->at(ss), fields, &ZN, &EM);
        calculateF(params, &IONS->at(ss), &ZN, &EM, &F);
        dZ4 = F*DT;

        // Assemble RK4 solution:
        Z1 = Z0 + (dZ1 + 2*dZ2 + 2*dZ3 + dZ4)/6;

				if ( isnan(ZN(0)) || isnan(ZN(1)) || isnan(ZN(2)) )
				{
					cout << "Z0(0) = " << Z0(0) << endl;
					cout << "Z1(0) = " << Z1(0) << endl;
					cout << "Z1(1) = " << Z1(1) << endl;
					cout << "Z1(2) = " << Z1(2) << endl;
				}

				// Interpolate fields at new particle position:
				interpEM(params, mesh, &IONS->at(ss), fields, &Z1, &EM);

				// Assign solution to output vector:
				switch (params->advanceParticleMethod)
				{
				case 1:
					// Do nothing, RK4 solved for (x, vpar, vper)
					break;
				case 2:
					// Calculate vper since RK4 solved for (x, vpar, mu)
					B  = EM(1);
					double mu = Z1(2);
					double vper = sqrt(2*mu*B/Ma);
					Z1(2)  = vper;
					break;
				}
        // End of RK solution:
        //==============================================================

        // Update new particle states:
        IONS->at(ss).x_p(ii)   = Z1(0);
        IONS->at(ss).v_p(ii,0) = Z1(1); // vpar
        IONS->at(ss).v_p(ii,1) = Z1(2); // vper
				IONS->at(ss).mu_p(ii)  = 0.5*Ma*pow(Z1(2),2)/EM(1) ; // mu

			} // End of parallel region
		}
	}//structure to iterate over all the ion species.
}

void PIC_TYP::assignCell(const params_TYP * params, const mesh_TYP * mesh, ions_TYP * IONS)
{
	// Total number of computational particles:
	int N_CP(IONS->N_CP_MPI);

	// Clear assignment function:
    IONS->wxc.zeros();
    IONS->wxl.zeros();
    IONS->wxr.zeros();

	#pragma omp parallel for default(none) shared(mesh, IONS, params, std::cout) firstprivate(N_CP)
    for(int ii=0; ii<N_CP; ii++)
    {
		// Calculate nearest grid point:
		double x_p     = IONS->x_p(ii);
		double x_p_min = params->mesh_params.Lx_min;
		double dx      = params->mesh_params.dx;
		signed int m = round( 0.5 + (x_p - x_p_min)/dx ) - 1;

		// Correct "m" near boundaries if out of bound:
		if ( m >= params->mesh_params.Nx)
		{
			cout << "Main AssignCell out of bound, m = " << m << endl;
			m = params->mesh_params.Nx - 1;
		}
		if ( m < 0)
		{
			cout << "Main AssignCell out of bound, m = " << m << endl;
			m = 0;
		}

		// Assign nearest grid point:
		IONS->mn(ii) = m;

		// Distance to nearest grid point:
		double x = mesh->xm(m) - x_p;

		// Assignment function:
		IONS->wxl(ii) = 0.5*pow(1.5 + ((x - dx)/dx),2); // Left:
		IONS->wxc(ii) = 0.75 - pow(x/dx,2);             // Center:
		IONS->wxr(ii) = 0.5*pow(1.5 - ((x + dx)/dx),2); // Right:

	} // parallel omp

}

void PIC_TYP::assignCell_AllSpecies(const params_TYP * params, const mesh_TYP * mesh, vector<ions_TYP> * IONS)
{
	// Iterate over all ion species:
  // =============================
  for(int ss=0;ss<IONS->size();ss++)
  {
    // Assign cell and calculate partial ion moments:
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
      // Assign cell:
			PIC_TYP::assignCell(params, mesh, &IONS->at(ss));
    }
	}
}

void PIC_TYP::extrapolateMoments_AllSpecies(const params_TYP * params, fields_TYP * fields, vector<ions_TYP> * IONS)
{
	// Iterate over all ion species:
  // =============================
  for(int ss=0;ss<IONS->size();ss++)
  {
  	// Assign cell and calculate partial ion moments:
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
			//Calculate partial moments:
			calculateIonMoments(params,fields,&IONS->at(ss));

			// Reduce IONS moments to PARTICLE ROOT:
			// =====================================
			MPI_ReduceVec(params, &IONS->at(ss).n_m);
			MPI_ReduceVec(params, &IONS->at(ss).nv_m);
			MPI_ReduceVec(params, &IONS->at(ss).P11_m);
			MPI_ReduceVec(params, &IONS->at(ss).P22_m);

			// Broadcast to all PARTICLE ranks:
			// ================================
			MPI_Bcast(IONS->at(ss).n_m.memptr()  , IONS->at(ss).n_m.size()  , MPI_DOUBLE, 0, params->mpi.COMM);
			MPI_Bcast(IONS->at(ss).nv_m.memptr() , IONS->at(ss).nv_m.size() , MPI_DOUBLE, 0, params->mpi.COMM);
			MPI_Bcast(IONS->at(ss).P11_m.memptr(), IONS->at(ss).P11_m.size(), MPI_DOUBLE, 0, params->mpi.COMM);
			MPI_Bcast(IONS->at(ss).P22_m.memptr(), IONS->at(ss).P22_m.size(), MPI_DOUBLE, 0, params->mpi.COMM);

			// Apply smoothing:
			// ===============
			for (int jj=0; jj<params->filtersPerIterationIons; jj++)
			{
				  smooth(&IONS->at(ss).n_m  , params->smoothingParameter);
				  smooth(&IONS->at(ss).nv_m , params->smoothingParameter);
				  smooth(&IONS->at(ss).P11_m, params->smoothingParameter);
				  smooth(&IONS->at(ss).P22_m, params->smoothingParameter);
				}

			// Calculate derived ion moments: Tpar_m, Tper_m:
			// ==============================================
			calculateDerivedIonMoments(params, &IONS->at(ss));
    }

		// 0th moment at various time levels are sent to fields processes:
    // =============================================================
		// Ion density:
    MPI_SendVec(params, &IONS->at(ss).n_m);
    MPI_SendVec(params, &IONS->at(ss).n_m_);
    MPI_SendVec(params, &IONS->at(ss).n_m__);
    MPI_SendVec(params, &IONS->at(ss).n_m___);
	}
}

void PIC_TYP::calculateIonMoments(const params_TYP * params, fields_TYP * fields, ions_TYP * IONS)
{
	// Ion density:
	IONS->n_m___ = IONS->n_m__;
	IONS->n_m__  = IONS->n_m_;
	IONS->n_m_   = IONS->n_m;

	// Ion flux:
	IONS->nv_m__ = IONS->nv_m_;
	IONS->nv_m_  = IONS->nv_m;

	// Calculate ion moments:
	eim(params,fields,IONS);
}

void PIC_TYP::eim(const params_TYP * params, fields_TYP * fields, ions_TYP * IONS)
{
	// Number of particles:
	int N_CP(IONS->N_CP_MPI);

	// Ion mass:
	double Ma = IONS->M;

	// Reference magnetic field:
	double B0 = params->mesh_params.B0; // Maybe use the current value at the reference location

	// Clearing content of ion moments:
	// ===============================
	IONS->ncp_m.zeros();
	IONS->n_m.zeros();
	IONS->nv_m.zeros();
	IONS->P11_m.zeros();
	IONS->P22_m.zeros();

	#pragma omp parallel default(none) shared(params, IONS, B0, Ma) firstprivate(N_CP)
	{
		// Create private moments:
		// ======================
		arma::vec ncp = zeros(params->mesh_params.Nx + 4);
		arma::vec n   = zeros(params->mesh_params.Nx + 4);
		arma::vec nv  = zeros(params->mesh_params.Nx + 4);
		arma::vec P11 = zeros(params->mesh_params.Nx + 4);
		arma::vec P22 = zeros(params->mesh_params.Nx + 4);

		// Assemble moments:
		// =================
		#pragma omp for
		for(int ii=0; ii<N_CP; ii++)
		{
			// Nearest grid point:
			int ix = IONS->mn(ii) + 2;

			// Particle velocity:
			double vpar = IONS->v_p(ii,0);
			double vper = IONS->v_p(ii,1);

			// vx component:
			arma::vec phi = 2*M_PI*randu<vec>(1);
			double vy = vper*cos(phi(0));

			// Particle weight:
			double a = IONS->a_p(ii);

			// Computational particle density per MPI:
			ncp(ix-1) += IONS->wxl(ii);
			ncp(ix)   += IONS->wxc(ii);
			ncp(ix+1) += IONS->wxr(ii);

			// Density:
			n(ix-1) += IONS->wxl(ii)*a;
			n(ix)   += IONS->wxc(ii)*a;
			n(ix+1) += IONS->wxr(ii)*a;

			// Particle flux density:
			nv(ix-1) += IONS->wxl(ii)*a*vpar;
			nv(ix) 	 += IONS->wxc(ii)*a*vpar;
			nv(ix+1) += IONS->wxr(ii)*a*vpar;

			// Stress tensor P11:
			P11(ix-1) += IONS->wxl(ii)*a*Ma*pow(vpar,2);
			P11(ix)   += IONS->wxc(ii)*a*Ma*pow(vpar,2);
			P11(ix+1) += IONS->wxr(ii)*a*Ma*pow(vpar,2);

			// Stress tensor P22:
			P22(ix-1) += IONS->wxl(ii)*a*Ma*pow(vy,2);
			P22(ix)   += IONS->wxc(ii)*a*Ma*pow(vy,2);
			P22(ix+1) += IONS->wxr(ii)*a*Ma*pow(vy,2);
		}

		// Reduce partial moments from each thread:
		// ========================================
		#pragma omp critical (update_ion_moments)
		{
			IONS->ncp_m.subvec(1,params->mesh_params.Nx) += ncp.subvec(2,params->mesh_params.Nx + 1);
			IONS->n_m.subvec(1,params->mesh_params.Nx)   += n.subvec(2,params->mesh_params.Nx + 1);
			IONS->nv_m.subvec(1,params->mesh_params.Nx)  += nv.subvec(2,params->mesh_params.Nx + 1);
			IONS->P11_m.subvec(1,params->mesh_params.Nx) += P11.subvec(2,params->mesh_params.Nx + 1);
			IONS->P22_m.subvec(1,params->mesh_params.Nx) += P22.subvec(2,params->mesh_params.Nx + 1);
		}

	}//End of the parallel region

	// Ghost contributions:
	// ====================
	fill4Ghosts(&IONS->ncp_m);
	fill4Ghosts(&IONS->n_m);
	fill4Ghosts(&IONS->nv_m);
	fill4Ghosts(&IONS->P11_m);
	fill4Ghosts(&IONS->P22_m);

	// Apply compression factor:
	// ========================
	arma::vec c = fields->Bx_m/B0;
	IONS->n_m   = IONS->n_m%c;
	IONS->nv_m  = IONS->nv_m%c;
	IONS->P11_m = IONS->P11_m%c;
	IONS->P22_m = IONS->P22_m%c;

	// Scale:
	// =====
	double A  = params->mesh_params.A0;
	double dx = params->mesh_params.dx;
	IONS->ncp_m /= params->mesh_params.dx;
	IONS->n_m   *= (1/A)*IONS->K/dx;
	IONS->nv_m  *= (1/A)*IONS->K/dx;
	IONS->P11_m *= (1/A)*IONS->K/dx;
	IONS->P22_m *= (1/A)*IONS->K/dx;
}

void PIC_TYP::calculateDerivedIonMoments(const params_TYP * params, ions_TYP * IONS)
{
	double Ma(IONS->M);

	// Ion pressures:
	arma::vec Ppar = IONS->P11_m - (Ma*IONS->nv_m % IONS->nv_m/IONS->n_m);
	arma::vec Pper = IONS->P22_m;

	// Ion temperatures:
	IONS->Tpar_m = Ppar/(F_E_DS*IONS->n_m);
	IONS->Tper_m = Pper/(F_E_DS*IONS->n_m);
}
