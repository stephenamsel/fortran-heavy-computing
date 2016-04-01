      subroutine mechanisms_run(omatrix, grid)
	  
     
	  !Produce operator-matrix of the following format:
	  !Index 1: blockID
	  !Index 2: Variable
	  !Index 3: Position in block of affected cell
	  !Index 4: Index in offset-list of offset of acting cell relative to affected cell
	  
	  !Inter-variable interactions are handled assuming that the acting variable is constant throughout the time-step
	  !These should be kept to a minimum 
	  
	  !Clear omatrix
	  
	  !Run model-equations here
      call pressure(grid, omatrix)
      call cooling(grid, omatrix)
      call heat_diffusion(grid, omatrix)
	  
	!Adjust for non-Cartesian coordinate-system here
	
      call geometry(omatrix)
	  
      end subroutine mechanisms_run
	  
      subroutine mechanisms_init()
	  !Sets parameters of omatrix, offset_list, actvar
	  !Parameters of omatrix determine its size
	  !Parameters of offset_list determine the locations of the acting elements of the matrix
	  
	  !Format of offset_list(nvar, 2, j):
	  !Index 1: variable affected
	  !Index 2, 3: {acting position, acting variable}
	  
      use mechanism_data
	  
      implicit none
	  
      integer mech, nmechs, noff, norepeat
      integer i, j, add, norepvar, offnow
      integer off, var2, istart, iend
	  
      integer, allocatable :: offset_raw(:,:,:)
      integer, allocatable :: raw_to_list(:)
	  
	  !Open config-file to get mechanism-data
      character(len = 200) line
	  !No more than 200 characters permitted on a line in config file
      write(len) "not started yet"
	  
      open(20, FILE = 'model_config')
	     do while (.true.) 
		    read(20, *) line
			
			!Exit at end of file
		    if line(1:4) == "Done" exit
			
			!Get size for offset_list
		    if line(1:6) == "offset" then
			
			   noff += 1
		    end if
	     end do
		 
	     rewind(20)
		 
	     allocate(offset_raw(nvar, 2, noff)) 
	     allocate(raw_to_list(noff)) 
		 
		 !Now get offsets into offset_raw
		 
	     offnow = 1
	     do while (.true.) 
		 
		    if line(1:6) == "offset" then
			    offactor = 1
			    do i=1, ndim
				   istart = 8 + (i-1) * 3
				   iend = 8 + (i) * 3 - 1
			       read(line(istart:iend), *) offdim
				   off += offdim * offactor
				   offactor *= blocksize(i)
			    end do
			   offset_raw(nvar, 1, offnow) = off
		    end if
			
		    if line(1:6) == "offvar" then
			    read(line(8:12), *) var2
			    do i=1, size(varlist)
				   if (varlist(i) == var2) offset_raw(nvar, 2, offnow)
			    end do
			    offnow += 1
		    end if

	     end do
	     close(20)
		 !Check for repetition, get mapping from raw to offset list
		 
	     norepeat = 0
	     do var = 1, nvar !FIX THIS!!!!
		    norepvar = 0 
	        do i = 1, noff
			   add = 1
			   do j = 1, i-1
				  if (offset_raw(var, 1, j) == offset_raw(var, 1, i) &
				     offset_raw(var, 1, j) == offset_raw(var, 1, i)) &
					 then
					 raw_to_list(i) = j
					 add = 0
				     exit
				  end if
			   end do
		       norepvar += add
		    end do
		    if (norepvar > norepeat) norepeat = norepvar
	     end do
		 
	     allocate(offset_list(nvar, 2, norepeat))
	  
      end subroutine mechanisms_init
	  
	  
      subroutine mechanisms_clear
	  !Frees data that had been allocated to run mechanisms and tore the operator matrix
	  
      end subroutine mechanisms_clear