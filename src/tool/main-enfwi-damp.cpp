extern "C" {
#include <rsf.h>
}

#ifdef _OPENMP
#include <omp.h>
#endif

#include <mpi.h>
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "logger.h"
#include "shot-position.h"
#include "damp4t10d.h"
#include "sf-velocity-reader.h"
#include "ricker-wavelet.h"
#include "essfwiframework.h"
#include "shotdata-reader.h"
#include "sfutil.h"
#include "sum.h"
#include "Matrix.h"
#include "updatevelop.h"
#include "updatesteplenop.h"
#include "enkfanalyze.h"
#include "common.h"
#include "environment.h"
#include "random-code.h"
#include "encoder.h"

namespace {
class Params {
public:
  Params();
  ~Params();
  void check();

public:
  sf_file vinit;        /* initial velocity model, unit=m/s */
  sf_file vreal;        /* real velocity model, unit=m/s */
  sf_file shots;        /* recorded shots from exact velocity model */
  sf_file vupdates;     /* updated velocity in iterations */
  sf_file absobjs;         /* absolute values of objective function in iterations */
  sf_file norobjs;         /* normalize values of objective function in iterations */
  int niter;            /* # of iterations */
  int nb;               /* size of the boundary */
  float vmin;
  float vmax;
  float maxdv;
  int nita;
  int seed;
  int nsample;
  int niterenkf;
  float sigfac;
  char *perin;

public: // parameters from input files
  int nz;
  int nx;
  float dz;
  float dx;
  int nt;
  int ng;
  int ns;
  float dt;
  float amp;
  float fm;
  int sxbeg;
  int szbeg;
  int gxbeg;
  int gzbeg;
  int jsx;
  int jsz;
  int jgx;
  int jgz;

public:
  int rank;
  int k;
  int np;
  int ntask; /// exactly the # of task each process owns
};

Params::Params() {
  vinit = sf_input ("vin");       /* initial velocity model, unit=m/s */
  vreal	= sf_input ("vreal");       /* real	velocity model, unit=m/s */
  shots = sf_input("shots");      /* recorded shots from exact velocity model */
  vupdates = sf_output("vout");   /* updated velocity in iterations */
  absobjs = sf_output("absobjs"); /* absolute values of objective function in iterations */
  norobjs = sf_output("norobjs"); /* normalized values of objective function in iterations */

  if (!sf_getint("niter", &niter)) { sf_error("no niter"); }      /* number of iterations */
  if (!sf_getfloat("maxdv", &maxdv)) sf_error("no maxdv");        /* max delta v update two iteration*/
  if (!sf_getint("nita", &nita))   { sf_error("no nita"); }       /* max iter refining alpha */
  if (!sf_getint("nsample", &nsample)){ sf_error("no nsample"); } /* # of samples for enkf */
  if (!sf_getint("niterenkf", &niterenkf)){ sf_error("no niterenkf"); } /* # of iteration between two enkf analyze */
  if (!(perin = sf_getstring("perin"))) { sf_error("no perin"); } /* perturbation file */
  if (!sf_getint("seed", &seed))   { seed = 10; }                 /* seed for random numbers */
  if (!sf_getfloat("sigfac", &sigfac))   { sf_error("no sigfac"); } /* sigma factor */

  /* get parameters from velocity model and recorded shots */
  if (!sf_histint(vinit, "n1", &nz)) { sf_error("no n1"); }       /* nz */
  if (!sf_histint(vinit, "n2", &nx)) { sf_error("no n2"); }       /* nx */
  if (!sf_histfloat(vinit, "d1", &dz)) { sf_error("no d1"); }     /* dz */
  if (!sf_histfloat(vinit, "d2", &dx)) { sf_error("no d2"); }     /* dx */
  if (!sf_histint(shots, "n1", &nt)) { sf_error("no nt"); }       /* total modeling time steps */
  if (!sf_histint(shots, "n2", &ng)) { sf_error("no ng"); }       /* total receivers in each shot */
  if (!sf_histint(shots, "n3", &ns)) { sf_error("no ns"); }       /* number of shots */
  if (!sf_histfloat(shots, "d1", &dt)) { sf_error("no dt"); }     /* time sampling interval */
  if (!sf_histfloat(shots, "amp", &amp)) { sf_error("no amp"); }  /* maximum amplitude of ricker */
  if (!sf_histfloat(shots, "fm", &fm)) { sf_error("no fm"); }     /* dominant freq of ricker */
  if (!sf_histint(shots, "sxbeg", &sxbeg)) {sf_error("no sxbeg");} /* x-begining index of sources, starting from 0 */
  if (!sf_histint(shots, "szbeg", &szbeg)) {sf_error("no szbeg");} /* x-begining index of sources, starting from 0 */
  if (!sf_histint(shots, "gxbeg", &gxbeg)) {sf_error("no gxbeg");} /* x-begining index of receivers, starting from 0 */
  if (!sf_histint(shots, "gzbeg", &gzbeg)) {sf_error("no gzbeg");} /* x-begining index of receivers, starting from 0 */
  if (!sf_histint(shots, "jsx", &jsx)) {sf_error("no jsx"); }       /* source x-axis  jump interval  */
  if (!sf_histint(shots, "jsz", &jsz)) { sf_error("no jsz"); }      /* source z-axis jump interval  */
  if (!sf_histint(shots, "jgx", &jgx)) { sf_error("no jgx"); }      /* receiver x-axis jump interval  */
  if (!sf_histint(shots, "jgz", &jgz)) { sf_error("no jgz"); }      /* receiver z-axis jump interval  */
  if (!sf_histint(shots,  "nb",&nb))        { sf_error("no nb"); }  /* thickness of sponge ABC  */
  if (!sf_histfloat(shots, "vmin", &vmin)) { sf_error("no vmin"); } /* minimal velocity in real model*/
  if (!sf_histfloat(shots, "vmax", &vmax)) { sf_error("no vmax"); } /* maximal velocity in real model*/
  //vmin -= 50;
  //vmax += 100;

  /**
   * output parameters
   */
  sf_putint(vupdates,   "n1", nz);
  sf_putint(vupdates,   "n2", nx);
  sf_putint(vupdates,   "n3", niter);
  sf_putfloat(vupdates, "d1", dz);
  sf_putfloat(vupdates, "d2", dx);
  sf_putint(vupdates,   "d3", 1);
  sf_putint(vupdates,   "o1", 0);
  sf_putint(vupdates,   "o2", 0);
  sf_putint(vupdates,   "o3", 0);
  sf_putstring(vupdates, "label1", "Depth");
  sf_putstring(vupdates, "label2", "Distance");
  sf_putstring(vupdates, "label3", "Iteration");
  sf_putint(absobjs, "n1", niter + 1);
  sf_putfloat(absobjs, "d1", 1);
  sf_putfloat(absobjs, "o1", 1);
  sf_putstring(absobjs, "label1", "Absolute");
  sf_putint(norobjs, "n1", niter + 1);
  sf_putfloat(norobjs, "d1", 1);
  sf_putfloat(norobjs, "o1", 1);
  sf_putstring(norobjs, "label1", "Normalize");

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  k = std::ceil(nsample * 1.0 / np);
  ntask = std::min(k, nsample - rank*k);

  check();
}

Params::~Params() {
  sf_close();
}

void Params::check() {
  if (!(sxbeg >= 0 && szbeg >= 0 && sxbeg + (ns - 1)*jsx < nx && szbeg + (ns - 1)*jsz < nz)) {
    sf_warning("sources exceeds the computing zone!\n");
    exit(1);
  }

  if (!(gxbeg >= 0 && gzbeg >= 0 && gxbeg + (ng - 1)*jgx < nx && gzbeg + (ng - 1)*jgz < nz)) {
    sf_warning("geophones exceeds the computing zone!\n");
    exit(1);
  }
}


std::vector<Velocity *> createVelDB(const Velocity &vel, const char *perin, int N, float dx, float dt) {
  std::vector<Velocity *> veldb(N);   /// here is all the velocity samples resides, others are pointer to this

  int modelSize = vel.nx * vel.nz;

  TRACE() << "add perturbation to initial velocity";
  std::ifstream ifs(perin);
  if (!ifs) {
    ERROR() << "cannot open file: " << perin;
    exit(EXIT_FAILURE);
  }

  std::vector<float> tmp(modelSize);
  std::vector<float> velOrig = vel.dat;
  std::transform(velOrig.begin(), velOrig.end(), velOrig.begin(), boost::bind(velRecover<float>, _1, dx, dt));

  for (int iv = 0; iv < N; iv++) {
    std::vector<float> ret(modelSize);
    ifs.read(reinterpret_cast<char *>(&tmp[0]), modelSize * sizeof(tmp[0]));
    std::transform(tmp.begin(), tmp.end(), velOrig.begin(), ret.begin(), std::plus<float>());
    std::transform(ret.begin(), ret.end(), ret.begin(), boost::bind(velTrans<float>, _1, dx, dt));
    veldb[iv] = new Velocity(ret, vel.nx, vel.nz);
  }

  ifs.close();

  return veldb;
}

std::vector<Velocity *> pCreateVelDB(const Velocity &vel, const char *perin, int N, float dx, float dt, int rank, int size) {
  std::vector<Velocity *> veldb(N);   /// here is all the velocity samples resides, others are pointer to this

  int modelSize = vel.nx * vel.nz;
	MPI_File fh;
	MPI_Offset file_size, offset;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_File_open(MPI_COMM_WORLD, perin, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Status status;
	offset = modelSize * sizeof(float) * N * rank;

  TRACE() << "parallel: add perturbation to initial velocity";
  if (!fh) {
    ERROR() << "cannot open file: " << perin;
    exit(EXIT_FAILURE);
  }

  std::vector<float> tmp(modelSize);
  std::vector<float> velOrig = vel.dat;
  std::transform(velOrig.begin(), velOrig.end(), velOrig.begin(), boost::bind(velRecover<float>, _1, dx, dt));

	int readSize = modelSize;

	printf("%d %d %d\n", modelSize, readSize, offset);

  for (int iv = 0; iv < N; iv++) {
    std::vector<float> ret(modelSize);
		MPI_File_seek(fh, offset, MPI_SEEK_SET);
		offset += readSize * sizeof(float);
		MPI_File_read(fh, &tmp[0], readSize, MPI_FLOAT, &status);
    std::transform(tmp.begin(), tmp.end(), velOrig.begin(), ret.begin(), std::plus<float>());
    std::transform(ret.begin(), ret.end(), ret.begin(), boost::bind(velTrans<float>, _1, dx, dt));
    veldb[iv] = new Velocity(ret, vel.nx, vel.nz);

  }

	MPI_File_close(&fh);

  return veldb;
}

std::vector<float *> generateVelSet(std::vector<Velocity *> &veldb) {
  std::vector<float *> velSet(veldb.size());
  for (size_t iv = 0; iv < veldb.size(); iv++) {
    velSet[iv] = &veldb[iv]->dat[0];
  }

  return velSet;
}

float calobj(const Damp4t10d &fmMethod, const std::vector<float> wlt,
    std::vector<float> &dobs, int ns, int ng, int nt) {
  // create random codes
	int seed = 1;
	RandomCodes r(seed);
  const std::vector<int> encodes = r.genPlus1Minus1(ns);

  Encoder encoder(encodes);
  std::vector<float> encsrc  = encoder.encodeSource(wlt);
  std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);

  std::vector<float> dcal(nt * ng, 0);
  fmMethod.EssForwardModeling(encsrc, dcal);
  fmMethod.removeDirectArrival(&encobs[0]);
  fmMethod.removeDirectArrival(&dcal[0]);

  std::vector<float> vsrc(nt * ng, 0);
  vectorMinus(encobs, dcal, vsrc);
  float obj = cal_objective(&vsrc[0], vsrc.size());

  return obj;
}

void scatterVelocity(std::vector<Velocity *> &veldb, const std::vector<Velocity *> &totalveldb, const Params &params) {
  int N = params.nsample;
  int rank = params.rank;
  int k = params.k;
  int ntask = params.ntask;
  int modelsize = veldb[0]->dat.size();

  if (rank == 0) { /// sender

    /// say N=7, np=3, k=ceil(N/np)=3, then each process owns the following velocity
    /// rank 0: 0, 1, 2   (tag: 0, 1, 2)
    /// rank 1: 3, 4, 5   (tag: 0, 1, 2)
    /// rank 2: 6         (tag: 0)
    /// each process except the last rank should own k samples. the last rank should own
    /// (N-(rank)*k) samples. (rank == np-1 for the last rank)
    /// the root process will send velocity tagged from 0 to k
    ///
    /// first prepare the velocity for the root process itself
    for (int i = 0; i < k; i++) {
      std::copy(totalveldb[i]->dat.begin(), totalveldb[i]->dat.end(), veldb[i]->dat.begin());
      DEBUG() << format("rank %d, vel%d %.20f") % rank % i % sum(veldb[i]->dat);
    }

    /// send velocity for other process, so sample is counting from k
    for (int isample = k; isample < N; isample++) {
      int rcvrank = isample / k;
      int tag = isample % k;
      MPI_Send(&totalveldb[isample]->dat[0], modelsize, MPI_FLOAT, rcvrank, tag, MPI_COMM_WORLD);
    }

  } else {
    /// (N-rank*k) for last rank process
    /// k for other process
    int sendrank = 0;
    for (int tag = 0; tag < ntask; tag++) {
      MPI_Status status;
      MPI_Recv(&veldb[tag]->dat[0], modelsize, MPI_FLOAT, sendrank, tag, MPI_COMM_WORLD, &status);

      DEBUG() << format("rank %d recv from rank %d, tag %d,  sum %.20f") % rank % sendrank % tag % sum(veldb[tag]->dat);
    }
  } //// end of dispatching velocity
}

void gatherVelocity(std::vector<Velocity *> &totalveldb, const std::vector<Velocity *> &veldb, const Params &params) {
  int N = params.nsample;
  int rank = params.rank;
  int k = params.k;
  int ntask = params.ntask;
  int modelsize = veldb[0]->dat.size();

  /// collect all the data from other process to rank 0
   if (rank != 0) {
     for (int tag = 0; tag < ntask; tag++) {
       DEBUG() << "1 rank " << rank << " tag: " << tag;
       MPI_Send(&veldb[tag]->dat[0], modelsize, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
       DEBUG() << format("2 send rank %d, tag %d, sum %.20f") % rank % tag % sum(veldb[tag]->dat);
     }
   } else { /// rank == 0
     /// collect from rank 0 itself first
     for (int i = 0; i < k; i++) {
       totalveldb[i]->dat = veldb[i]->dat;
     }

     // receive from other process
     // send velocity for other process, so sample is counting from k
     for (int isample = k; isample < N; isample++) {
       int sendrank = isample / k;
       int tag = isample % k;
       DEBUG() << format("isample %d") % isample;
       DEBUG() << totalveldb[isample]->dat[0];

       MPI_Status status;
       MPI_Recv(&totalveldb[isample]->dat[0], modelsize, MPI_FLOAT, sendrank, tag, MPI_COMM_WORLD, &status);
       DEBUG() << format("3 recv from rank %d, tag %d, sum %.20f") % sendrank % tag % sum(totalveldb[isample]->dat);
     }

     for (size_t i = 0; i < totalveldb.size(); i++) {
       DEBUG() << format("after recv correct vel%d %.20f") % i % sum(totalveldb[i]->dat);
     }

   }
}

void slownessL1L2Norm(const float *accurate2, const float *curr, const Params &params, float &l1norm, float &l2norm) {
	int bz = params.nb;
	int bx = bz;
	int nx = params.nx;
	int nz = params.nz;

	float norm1 = 0;
	float norm2 = 0;
	float sumSquare = 0;
	float sum = 0;

	std::vector<float> accurate(accurate2, accurate2 + nx * nz);
	std::transform(accurate.begin(), accurate.end(), accurate.begin(), boost::bind(velRecover<float>, _1, params.dx, params.dt));
	std::vector<float> tmp(curr, curr + nx * nz);
	std::transform(tmp.begin(), tmp.end(), tmp.begin(), boost::bind(velRecover<float>, _1, params.dx, params.dt));

	for (int ix = bx; ix < nx - bx; ix++) {
		for (int iz = 0; iz < nz - bz; iz++) {
			int offset = ix * nz + iz;

			float accvel = 1.0f / accurate[offset]; /// slowness
			float curvel = 1.0f / tmp[offset];      /// slowness
			float absdiff = std::abs(accvel - curvel);

			norm2     += absdiff * absdiff;
			sumSquare += accvel * accvel;

			norm1  += absdiff;
			sum    += accvel;
		}
	}

	l1norm = norm1 / sum;
	l2norm = norm2 / sumSquare;
}
}/// end of name space


int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  sf_init(argc, argv); /* initialize Madagascar */
  Environment::setDatapath();

  Params params;
  /// configure logger
  std::string logfile = std::string("enfwi-damp-") + boost::lexical_cast<std::string>(params.rank) + ".log";
  FILELog::setLogFile(logfile);

  int nz = params.nz;
  int nx = params.nx;
  int nb = params.nb;
  int ng = params.ng;
  int nt = params.nt;
  int ns = params.ns;
  float dt = params.dt;
  float fm = params.fm;
  float dx = params.dx;
  float vmin = params.vmin;
  float vmax = params.vmax;
  float maxdv = params.maxdv;
  float sigfac = params.sigfac;
  int nita = params.nita;
  int N = params.nsample;
  int niterenkf = params.niterenkf;
  int rank = params.rank;
  int k = params.k;
  int ntask = params.ntask;

  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  ShotPosition allGeoPos(params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);
  Damp4t10d fmMethod(allSrcPos, allGeoPos, dt, dx, fm, nb, nt);
  std::vector<float> wlt = rickerWavelet(nt, fm, dt, params.amp);

  /// read and broadcast velocity
  SfVelocityReader velReader(params.vinit);
  Velocity v0(nx, nz);
  velReader.readAndBcast(&v0.dat[0], nx * nz, 0);
  Velocity exvel = fmMethod.expandDomain(v0);
  fmMethod.bindVelocity(exvel);

  SfVelocityReader velReader_real(params.vreal);
  Velocity v0_real(nx, nz);
  velReader_real.readAndBcast(&v0_real.dat[0], nx * nz, 0);
  Velocity exvel_real = fmMethod.expandDomain(v0_real);

  std::vector<float> absobj;
  std::vector<float> norobj;

  /// read and broadcast dobs
  std::vector<float> dobs(ns * nt * ng);     /* all observed data */
  if (rank == 0) {
    ShotDataReader::serialRead(params.shots, &dobs[0], ns, nt, ng);
  }
  MPI_Bcast(&dobs[0], dobs.size(), MPI_FLOAT, 0, MPI_COMM_WORLD);

  UpdateVelOp updatevelop(vmin, vmax, dx, dt);

  //std::vector<Velocity *> totalveldb;  /// only for rank 0
  std::vector<Velocity *> veldb(ntask); /// each process owns # of velocity


	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int model_size = veldb[0]->nx * veldb[0]->nz;
  //std::vector<Velocity *> veldb2(ntask); /// each process owns # of velocity
	veldb = pCreateVelDB(exvel, params.perin, veldb.size(), dx, dt, rank, size);
  //EnkfAnalyze enkfAnly2(fmMethod, wlt, dobs, sigfac);

	/*
  if (rank == 0) { /// sender
    totalveldb = createVelDB(exvel, params.perin, N, dx, dt);
    for (size_t i = 0; i < totalveldb.size(); i++) {
      DEBUG() << format("correct vel%d %.20f") % i % sum(totalveldb[i]->dat);
    }
    DEBUG();
  }
  for (size_t iv = 0; iv < veldb.size(); iv++) {
    veldb[iv] = new Velocity(exvel.nx, exvel.nz);
  }
	*/

  /// it is a bad implementation. one process send each velocity to all other processes
  /// if the number of samples become large, it would be the bottleneck.
  /// TODO: refactor the Velocity class. keep the real data of all velocity in ONE array and
  /// each Velocity has a pointer that points to the correct location.

  //TODO: need modifying, delete scatter and gather, make processes reading files in parallel
  //scatterVelocity(veldb, totalveldb, params);
  DEBUG() << "dispatching velocity finished!";
//  MPI_Barrier(MPI_COMM_WORLD);

	/*
	for(int i = 0 ; i < veldb.size() ; i ++)
		enkfAnly2.check(veldb[i]->dat, veldb2[i]->dat);
		*/
//	MPI_Barrier(MPI_COMM_WORLD);


  std::vector<Damp4t10d *> fms(ntask);
  std::vector<UpdateSteplenOp *> usl(ntask);
  std::vector<EssFwiFramework *> essfwis(ntask);

  for (size_t i = 0; i < essfwis.size(); i++) {
    fms[i] = new Damp4t10d(fmMethod);
    fms[i]->bindVelocity(*veldb[i]);

    usl[i] = new UpdateSteplenOp(*fms[i], updatevelop, nita, maxdv);
    essfwis[i] = new EssFwiFramework(*fms[i], *usl[i], updatevelop, wlt, dobs);
  }

  EnkfAnalyze enkfAnly(fmMethod, wlt, dobs, sigfac);


  /// collect all the data from other process to rank 0
  //gatherVelocity(totalveldb, veldb, params);

  std::vector<float *> velset = generateVelSet(veldb);
  //std::vector<float *> totalVelSet;

	/*
  if (rank == 0) {
    totalVelSet = generateVelSet(totalveldb);
  }
	*/

  //enkfAnly.analyze(totalVelSet, velset);

  /// for regulalization
  float initLambdaRatio = 0.5;
  float local_n = velset.size();
  Matrix ratioSet(local_n, 2);  /// 0 for muX, 1 for muZ
  Matrix lambdaSet(local_n, 2); /// 0 for lambdaX, 1 for lambdaZ
  std::fill(ratioSet.getData(), ratioSet.getData() + ratioSet.size(), initLambdaRatio);
  enkfAnly.initLambdaSet(velset, lambdaSet, ratioSet);
  enkfAnly.pAnalyze(velset, lambdaSet, ratioSet);


//  enkfAnly.pAnalyze(velset);


  //TODO: need modifying, createAMean
	std::vector<float> vvt = enkfAnly.pCreateAMean(velset, N);

  if (rank == 0) {
    /// calculate objective function

    //std::vector<float> vv = enkfAnly.createAMean(totalVelSet);
    std::vector<float> vv = vvt;
		//enkfAnly.check(vvt, vv);
		//exit(1);
    Velocity newvel(vv, fmMethod.getnx(), fmMethod.getnz());
    fmMethod.bindVelocity(newvel);
    float obj = calobj(fmMethod, wlt, dobs, ns, ng, nt);
    absobj.push_back(obj);
    norobj.push_back(obj / absobj[0]);
  }

  /// after enkf, we should scatter velocities
  //scatterVelocity(veldb, totalveldb, params);




  srand(params.seed + params.rank);
  TRACE() << "iterate the remaining iteration";
  for (int iter = 0; iter < params.niter; iter++) {
    TRACE() << "FWI for each velocity";
    DEBUG() << "\n\n\n\n\n\n\n";

    for (size_t ivel = 0; ivel < essfwis.size(); ivel++) {
      int absvel = rank * k + ivel + 1;
      INFO() << format("iter %d, rank %d on %dth velocity, sum %f") % iter % rank % absvel % sum(veldb[ivel]->dat);
      double lambdaX = lambdaSet.getData()[ivel * lambdaSet.getNumRow()];
      double lambdaZ = lambdaSet.getData()[ivel * lambdaSet.getNumRow() + 1];
      essfwis[ivel]->epoch(iter, lambdaX, lambdaZ);
    }

    TRACE() << "enkf analyze and update velocity";
    //gatherVelocity(totalveldb, veldb, params);
    if (iter % niterenkf == 0) {
			for(int ivel = 0 ; ivel < essfwis.size() ; ivel ++)
			{
				int absvel = rank * k + ivel + 1;
				INFO() << format("%4d/%d iter, velset[%d/%d] data rss: %g") % iter % params.niter % absvel % N % essfwis[ivel]->getUpdateObj();

				float l1norm, l2norm;
				slownessL1L2Norm(&exvel_real.dat[0], &veldb[ivel]->dat[0], params, l1norm, l2norm);
				INFO() << format("%4d/%d iter, velset[%d/%d] slowness l1norm: %g, slowness l2norm: %g") % iter % params.niter % absvel % N % l1norm % l2norm;
			}

      //enkfAnly.analyze(totalVelSet, velset);
			enkfAnly.pAnalyze(velset, lambdaSet, ratioSet);

    }

		std::vector<float> vvt = enkfAnly.pCreateAMean(velset, N);
    if (rank == 0) {
      /// output velocity
      //std::vector<float> vv = enkfAnly.createAMean(totalVelSet);
      std::vector<float> vv = vvt;

      fmMethod.sfWriteVel(vv, params.vupdates);

      /// calculate the objective function for the updated velocity
      Velocity newvel(vv, fmMethod.getnx(), fmMethod.getnz());
      fmMethod.bindVelocity(newvel);

      float obj = calobj(fmMethod, wlt, dobs, ns, ng, nt);
      INFO() << format("iter %d, objval for meam model %e") % iter % obj;

      absobj.push_back(obj);
      norobj.push_back(obj / absobj[0]);
    }
    //scatterVelocity(veldb, totalveldb, params);
  }

  /// write objective function values
  if (rank == 0) {
    sf_floatwrite(&absobj[0], absobj.size(), params.absobjs);
    sf_floatwrite(&norobj[0], norobj.size(), params.norobjs);
  }

  /// release memory
	/*
  if (rank == 0) {
    for (size_t i = 0; i < totalveldb.size(); i++) {
      delete totalveldb[i];
    }
  }
	*/
  for (int i = 0; i < ntask; i++) {
    delete veldb[i];
    delete fms[i];
    delete usl[i];
    delete essfwis[i];
  }


  MPI_Finalize();
  return 0;
}
