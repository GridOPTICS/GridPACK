#include "counters.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <stdint.h>

#include <mpi.h>
#include <tascel.h>
#include <tascel/UniformTaskCollectionSplit.h>

/****************Function that calls back into app*****************/

static int me;
extern "C" void updatecnt_(MPI_Fint *);

/* TODO: call the actual application function */
static void actual_func(int64_t taskid) {
//  printf("%d: Task id=%d\n", me, (int)taskid);

  MPI_Fint tid=taskid;
  updatecnt_(&tid);
  //do the actual work by calling the app routine
}


/***********External API called from app code***********************/

//forward declarations
using namespace tascel;
using std::vector;

static void amBarrier();
static void tascelInit(int intraRanks);
static void tascelFinalize();
static void task_func(tascel::UniformTaskCollection *utc, void *tt, int tt_len,
		      void *pldata, int pldata_len, vector<void *> data_bufs, int thd);
struct task_dscr_t {
  uint64_t taskid;
};
static ProcGroup* pgrp = NULL;
static UniformTaskCollectionSplit *utc = NULL;

void counters_create_() {
  tascelInit(1);
}

void counters_destroy_() {
  amBarrier();
  tascelFinalize();
}

void counters_process_(MPI_Fint *count) {
  MPI_Barrier(MPI_COMM_WORLD);
  /*leaving MPI world for tascel*/

  const int64_t nproc = theTwoSided().numProcRanks().toInt();
  me = theTwoSided().getProcRank().toInt();
  const int64_t total_ntasks = *count;
  assert(total_ntasks>=0);
  const int64_t ntasks_per_proc = (int64_t)ceil(1.0*total_ntasks/nproc); 
  
  const int64_t tasklo = ntasks_per_proc * me;
  const int64_t taskhi = min(tasklo+ntasks_per_proc,total_ntasks)-1;

#if 1
  {
    TslFuncRegTbl* frt  = new TslFuncRegTbl();
    TslFunc tf = frt->add(task_func);

    TaskCollProps props;
    props.functions(tf,frt)
      .taskSize(sizeof(task_dscr_t))
      .maxTasks(ntasks_per_proc*10);

    utc = new UniformTaskCollectionSplit(props, 0);

    amBarrier();
    //add the tasks

    task_dscr_t tdscr;
    for (int64_t i=tasklo; i<=taskhi; ++i) {
      tdscr.taskid = i;
      utc->addTask(&tdscr, sizeof(tdscr));
    }
    amBarrier();
    utc->process();
    amBarrier();
    utc->printStats();
    delete utc;
    utc = NULL;
    delete frt;
  }
  amBarrier();
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  /*going back to MPI world*/
}


/****************Internal routines that should not be affected********/

static void amBarrier() {
  int epoch = pgrp->signalBarrier();
  while(!pgrp->testBarrier(epoch)) {
    AmListenObjCodelet<NullMutex>* codelet;
    if((codelet=theAm().amListeners[0]->progress()) != NULL) {
      codelet->execute();
    }
  }
}

static void tascelInit(int intraRanks) {
  TascelConfig::initialize(intraRanks, MPI_COMM_WORLD);
  pgrp = ProcGroup::construct();
}

static void tascelFinalize() {
  delete pgrp;
  pgrp = NULL;
  TascelConfig::finalize();
}

static void task_func(UniformTaskCollection *utc, void *tt, int tt_len,
		      void *pldata, int pldata_len, vector<void *> data_bufs, int thd) {
  assert(tt != NULL);
  assert(utc!=NULL);
  assert(tt_len == sizeof(task_dscr_t));
  assert(pldata == NULL);
  assert(pldata_len == 0);
  assert(data_bufs.size()==0);
  assert(thd == 0);

  task_dscr_t& task = *(task_dscr_t*)tt;
  actual_func(task.taskid);
}

