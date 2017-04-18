module fem
  use typy  
  public :: solve_pde
  private :: terminal_info, exitme, close_all_observe
  integer(kind=ikind), private :: itcount, obs_pos, obs_count, nptimes
  logical, private :: go_clusters
  logical, private :: printtime


  contains

    !> solves the nonstationary pde problem in time
    subroutine solve_pde(success)
      use typy
      use globals
      use global_objs
      use pde_objs
      use femmat
      use feminittools
      use postpro
      use decomposer
      use decomp_tools
      use decomp_vars
      use core_tools
      use debug_tools

      
      logical, intent(out) :: success
      integer :: ierr
      integer :: file_itcounts
      integer(kind=ikind) :: i
      integer :: fileid
      real(4), dimension(3) :: act_time
      real(4) :: logtime=0.0, logtime_dt=900
      

      minimal_dt = dtmax
      time_step = init_dt
      if (.not. drutes_config%run_from_backup) time = 0
      itcum = 0
      printtime = .false.
      obs_count = ubound(observe_time,1)
      
      if (.not. drutes_config%run_from_backup) then
	obs_pos = 1
      else
	obs_pos = 0
	do i=1, ubound(observe_time,1)
	  if (observe_time(i)%value > time) then
	    obs_pos = i
	    EXIT
	  end if
	end do
	if (obs_pos == 0) then
	  obs_pos = ubound(observe_time,1) + 1
	end if
      end if


      call find_unit(file_itcounts)
      open(unit=file_itcounts, file="out/itcounts", action="write", status="replace")

      call find_unit(file_dt)
      open(unit=file_dt, file="out/dt", action="write", status="replace")
      
     call make_print("separately")
     
      call write_log("go 4 solving")
      
      do

        if (minimal_dt > time_step) then
          minimal_dt = time_step
        end if

	call pde_common%treat_pde(ierr, itcount, success)
	
	itcum = itcum + itcount

	if (success) then
	  write(unit=file_itcounts, fmt=*) time, itcount
	  call flush(file_itcounts)
	  write(unit=file_dt, fmt=*) time, time_step
	  call flush(file_dt)
	  

	  time = time + time_step

	  call write_obs()
	  if (printtime) then
	    do i=1, nptimes
	      	      
	      call make_print("separately", observe_time(obs_pos)%value, observe_time(obs_pos)%name)
	      
	      obs_pos = obs_pos + 1
	      write(unit=terminal, fmt=*)  " " //achar(27)//'[97m',"I: print time was reached, making output files..."&
					    //achar(27)//'[0m'
	
	    end do
	    if (drutes_config%it_method > 0) then
	      call print_domains("separately")
	      call print_elements_dd("separately")
	    end if
	    printtime = .false.
	  end if
	end if

	call evol_dt(ierr, itcount)

	call terminal_info(ierr)
	
	if (time >= end_time) then
	  call make_print()
	  call close_all_observe()
	  success = .true.
	  RETURN
	end if
	call etime(act_time(1:2), act_time(3))
	if (act_time(3) > logtime) then
	  call write_log("current simulation time:", real1=time, text2="total iteration count:", int2=itcum, hidden=.true.)
	  logtime = logtime + logtime_dt
	end if
	if (cpu_time_limit .and. act_time(3) >= cpu_max_time) then
	  call write_log("maximal allowed CPU time reached, exiting....")
	  success = .false.
	  RETURN
	end if
      end do

	  
    end subroutine solve_pde
    
    subroutine close_all_observe()
      use typy
      use pde_objs
      
      integer(kind=ikind) :: i, name

      do i=1, ubound(pde,1)
        do name=1, ubound(pde(i)%obspt_unit,1)
          close(pde(i)%obspt_unit(name))
        end do
      end do
    end subroutine close_all_observe
        

    subroutine evol_dt(ierr, itcount)
      use typy
      use globals
      use feminittools
      use core_tools
      use debug_tools
      use pde_objs
      
      integer, intent(in) :: ierr
      integer(kind=ikind), intent(in) :: itcount
      logical :: success_it


      select case(ierr)
	case(0)
	  dtprev = time_step
	  if (itcount < 0.25*max_itcount) then
	    time_step = min(dtmax, 1.06*time_step) 
	  else
	    time_step = time_step
	  end if
	  
	  if (itcount >0.45*max_itcount) then
	    time_step = max(dtmin, 0.9*time_step)
	  end if
	    
	  success_it = .true.
	case(1,2)
	  dtprev = time_step
	  time_step = max(dtmin, 0.8*time_step)
	  success_it = .false.
	  if (abs(time_step - dtmin) < epsilon(dtmin)) then
	    call exitme()
	  end if
      end select
      
  

      if (success_it) then	
	if (obs_pos <= obs_count) then
	  if (time + time_step >= observe_time(obs_pos)%value) then
	    select case(observe_info%method)
	      case(1)
		printtime = .true.
		nptimes = 1
		if (observe_time(obs_pos)%value - time > dtmin) then
		  time_step = observe_time(obs_pos)%value - time
		else
		  time_step = dtmin
		  call write_log(text=&
		    "Warning: Observation time adjusted due to print times, but the required time step seems to short:", &
		      real1=time_step)
		end if
		call write_log(text="Observation time adjusted due to print times, current time step is:", real1=time_step)
	      case(2)
		printtime = .true.
		nptimes = 1
		do 
		  if (obs_pos+nptimes > ubound(observe_time,1)) then
		    EXIT
		  else
		    if (time + time_step >= observe_time(obs_pos+nptimes)%value) then
		      nptimes = nptimes+1
		    else
		      EXIT
		    end if
		  end if
		end do
	    end select
	  end if
	end if
      end if

    end subroutine evol_dt

    subroutine terminal_info(ierr)
      use globals
      use printtools

      integer, intent(in) :: ierr

      write(unit=terminal, fmt=*) " "
      write(unit=terminal, fmt=*) " "
      if (.not. www) then
	write(unit=terminal, fmt=*) "current working directory is: " //achar(27)//'[101m', trim(dir_name) //achar(27)//'[0m'
      end if
      write(unit=terminal, fmt=*) " "


      select case(ierr)
	case(0)
	  write(unit=terminal, fmt=*)" " //achar(27)//'[93m',  "--------------------OK-------------------------------------" &
					  //achar(27)//'[0m'
	  write(unit=terminal, fmt=*) " " //achar(27)//'[92m', "actual simulation time: " //achar(27)//'[0m', time, "  | ", &
				      "previous time step: ", dtprev
	  write(unit=terminal, fmt=*) "       _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _"
          write(unit=terminal, fmt=*) "  "
	  write(unit=terminal, fmt=*) "total iteration count: ", itcum,"         | ", "current iteration count: ", itcount
	  write(unit=terminal, fmt=*) "  "
	  write(unit=terminal, fmt="(a, I4, a, I4, a)") " observation time files written: ", obs_pos-1, "/", obs_count
	  write(unit=terminal, fmt=*) "  "
	  call time2finish()
	  write(unit=terminal, fmt=*)" " //achar(27)//'[93m',  "-----------------------------------------------------------" &
					  //achar(27)//'[0m'
	case(1)
          write(unit=terminal, fmt=*)" " //achar(27)//'[91m',  "--------------------WARNING!-------------------------------" &
					  //achar(27)//'[0m'
          write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "slow convergence of the Picard method, time step decreased: " &
				      //achar(27)//'[0m', time_step

	  write(unit=terminal, fmt=*) "current simulation time:", time
	  write(unit=terminal, fmt=*) "total iteration count:", itcum
          write(unit=terminal, fmt=*)" " //achar(27)//'[91m',  "------------------------------------------------------------" &
					  //achar(27)//'[0m'
	case(2)
	  write(unit=terminal, fmt=*)" " //achar(27)//'[91m',  "--------------------INFO!-------------------------------" &
					  //achar(27)//'[0m'
          write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "RCZA method forced time step decrease: " &
				      //achar(27)//'[0m', time_step

	  write(unit=terminal, fmt=*) "current simulation time:", time
	  write(unit=terminal, fmt=*) "total iteration count:", itcum
          write(unit=terminal, fmt=*)" " //achar(27)//'[91m',  "------------------------------------------------------------" &
					  //achar(27)//'[0m'
      end select
    end subroutine terminal_info
    
    subroutine exitme()
      use typy
      use globals
      use feminittools
      use core_tools
      use postpro
      
      integer(kind=ikind) :: selection, i, ii, itold
      real(kind=rkind) :: valold
    
      call make_print("separately")
      if (www) then
	call write_log("failed to converge, you can run the code again with different setup, &
	the last output file can be used to relaunch your computation")
	ERROR STOP
      else
	call write_log("failed to converge, user can change some values now")
	print *, "select the following:"
	print *, "                     1 - update number of iterations"
	print *, "                     2 - update minimal time step"
	print *, "                     3 - update Picard iteration threshold"
	print *, "                     4 - exit the code"
! 	read(*,iostat=i) selection
	i = 0
	read(*,*) selection
	if (i == 0) then
	  select case(selection)
	    case(1)
	      print *, "the number of maximal number of iterations is:", max_itcount, "type new value now"
	      itold = max_itcount
	      do
! 		read(*, iostat=ii) max_itcount
		ii = 0
		read(*, *) max_itcount
		if (ii /= 0) then
		  print *, "that wasn't well done, try over again"
		else
		  exit
		end if
	      end do
	      
	      call write_log("the number of maximal number of iterations was updated for", int1=max_itcount,&
	      text2="the previous number of iterations was", int2 = itold)
	      
	      
	    case(2)
	      print *, "the minimal time step is:", dtmin, "type new value now"
	      valold = dtmin
	      do
! 		read(*, iostat=ii) dtmin
		ii = 0
		read(*, *) dtmin
		if (ii /= 0) then
		  print *, "that wasn't well done, try over again"
		else
		  exit
		end if
	      end do
				  
	      call write_log("the minimal time step was updated for", real1=dtmin,&
	      text2="the previous time step was", real2 = valold)
	    
	    case(3)
	      print *, "the Picard iteration criterion is:", iter_criterion, "type new value now"
	      valold = iter_criterion
	      do
! 		read(*, iostat=ii) iter_criterion
		ii = 0
		read(*, *) iter_criterion
		if (ii /= 0) then
		  print *, "that wasn't well done, try over again"
		else
		  exit
		end if
	      end do
	      call write_log("the iteration criterion for the Picard method was updated for", &
	      real1=iter_criterion,&
	      text2="the previous value was", real2 = valold)
	    case(4)
	    
	      call write_log("You have decided to EXIT this unstable and tormented computation. BYE!")
	      STOP
	    
	    case default
	      print *, "You have typed an unsupported keyword here, I'll see you here probably over again :)"
	      CONTINUE
	  end select
	else
	  print *, "You have typed an unsupported keyword here, I'll see you here probably over again :)"
	  CONTINUE
	end if
      end if
    end subroutine exitme
    
    



end module fem
