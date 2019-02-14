!* By Jia Shi 
!# control the procedures

program cg_evsl

  use mpi
  use omp_lib
  
  use para_mod,                 only: para_init,pin           
  use utility_mod    
  use geometry_mod 
  !
  use cg_models_mod,            only: cg_load_models
  use cg_create_matrix_mod,     only: cg_create_matrix

  use cg_matvec_mod 
  use pevsl_mod

  implicit none 

  integer                          :: ierr

  call mpi_init(IERR)

  call para_init("global_conf")
  call openlogfile(pin%comm)
  call report_time(pin%comm)
  call writelog("Starting main program",pin%comm)
  
  call writelog("Building geometry",pin%comm)
  call Build_geometry()
  call mpi_barrier(pin%comm,ierr)
  call report_time(pin%comm)
  
  call writelog("Loading models",pin%comm)
  call cg_load_models()
  call report_time(pin%comm)
  
  call writelog("Creating matrix",pin%comm)
  call cg_create_matrix()
  call report_time(pin%comm)

  call setupmatvec()
  call report_time(pin%comm)

  call pnm_apply_pevsl() 
  call writelog("solve the eigenvalue problem",pin%comm)
  call report_time(pin%comm)
  
  call writelog("finalize everything",pin%comm)
  call report_time(pin%comm)
  call MPI_FINALIZE(IERR)


end program cg_evsl 
  

