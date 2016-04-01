      subroutine explicit(dtexp, omatrix)
	  
      use grid_data
      use mechanism_data, only :: offset_list
		 
      implicit none
      integer lb, i, j, k
		 
		 
       do lb = 1, size(blocklist)
	     do i=1, product(blocksize)
		    do k = 1, nvars
			   var = varorder(k)
			   grid(blockID, var, i) = old_grid(blockID, var, i)
			   do j = 1, size(offset_list(var, :)) Can make it depend upon position if that matters
				  grid(BlockID, var, i) += newgrid(BlockID, var, &
				     i + offset_list(var, j)) * omatrix(blockID, var, i,j) * dtexp
				  
			   end do
		    end do
         end do
      end do
	  
	  
      end subroutine explicit
	  
      subroutine explicit()
	  
      use grid_data
      use mechanism_data :: nonzero_dist, nonzero
      
	  
      implicit none
	  
      integer :: i, j, k, n, nn, a, b
      logical :: skip
      integer, dimension(ndim) :: local_index, interaction_index
	  !Detect furthest non-zero element from the diagonal in each variable to limit number of operations
	  !MOVE THIS TO INITIALIZATION, USE IT TO LIMIT GRAPHS IN IMPLICIT SECTION?
	  
	  
	  !Multiply, limiting loops by what is furthest from the diagonal
      do i= 1, product(blocksize)
	     do n=1, ndim
		 !Convert to local indices
		    local_index(n) = mod(i, product(blocksize(n:ndim)))
	        if local_index(n) < block_boundary(1, n) .or. &
			   local_index(n) > block_boundary(2, n) then
			      skip = .true.
				  exit
		    end if
	     end do
	     if skip then
	        cycle
	     end if
	!Loop over variables for inter-variable interaction?
	     a = 1, size(offsets(var, :))
	
	     do a = 1,ndim
		    do b = 1,ndim
		    !Loop over elements within space defined by maximum distance of non-0 elements from diagonal
		
		!REPLACE NONZER-MATRIX WITH OFFSETS (INCLUDING OFFSETS FOR DIFFERENT DIMENSIONS)
		
		    interaction_index(:) = nonzero(a, b, 1,:) !lower bound of interacting space, element of variable b acting on a
			!ex: Pressure acting on velocity in hydrodynamics: b = P index, a = V_x index
		    do k = 1, product(nonzero_dist(a, b, :)) !nonzero_dist(a, b, i) = nonzero(2,i) - nonzero(1,i)
			   interaction_index(ndim) += 1
			   
			   !"nonzero" is in intuitive rectangular coordinates, not unwrapped vector. Get coordinate.
			   do nn = 1, ndim - 1
				  if interaction_index(ndim+1-nn) > nonzero(a, b, 2,ndim+1-nn) then
				     interaction_index(ndim + 1 - nn) = nonzero(a, b, 1,ndim+1-nn)
					 interaction_index(ndim - nn) += 1
				  end if
			   end do
			   !Convert to vector-coordinate
			   j = i
			   matrix_vector_index = 0
			   do nn = 1, ndim
			      j += interaction_index(nn) * &
					 product(blocksize(nn:ndim)) / blocksize(nn)
				  matrix_vector_index += (interaction_index(nn) - &
					 nonzero(a, b, 1,ndim+1-nn)) * &
					 product(blocksize(nn:ndim)) / blocksize(nn)
			   end do
			   !FOR OF OMATRIX: OMATRIX(BLOCKID, VAR, AFFECTED_INDEX, ACTING_INDEX)
			   !Do convolution-element
			   !!!!MODIFY THIS TO REDUCE MEMORY-REQUIREMENT!!!!
			   !CALL MECHANISMS HERE TO CALCULATE DYNAMICALLY?
			   grid(blockID, a, i) += 
			   omatrix(blockID, i, a, b, matrix_vector_index) * &
				  old_grid(blockID, j, b)
			   
		    end do
	     end do
			
	     end do
      end do
	  
      end subroutine explicit