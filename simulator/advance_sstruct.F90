      #include "variables.h"
	  !IMPLEMENT VARIABLES.H OR INCLUDE ITS DATA IN GRID_DATA
      #include "tree.h"
      #include <math.h>
      #include "_hypre_utilities.h"
      #include "HYPRE_sstruct_mv.h"
      #include "HYPRE_sstruct_ls.h"
      #include "HYPRE.h"
	  
      subroutine sstruct_init
      use sstruct_data
      use grid_data      
	  
      implicit none
      integer*8 :: hypre_grid !grid-structure information
      integer*8 :: hypre_graph !information defining non-zero part of operator matrix
      
      integer, dimension(2, ndim) :: extents
      integer, dimension(nvars) :: vartype
      integer, dimension(2 ** ndim, ndim, ndim) :: offcorner !Used for matching indices
      integer, dimension(2 ** ndim, ndim) :: corner
      integer, dimension(ndim, 2, 2, ndim) :: connector, off_connector
	  
	  !connector indices: dimension, low / high face, low / high index on face, coordinate
	  
      integer i, j, m
      integer nlref
	  
      call HYPRE_SStructFACCreate(MPI_Comm, solver)
      
      do i=1, ndim
         extents(1,i) = 1
         extents(2,i) = blocksize(i)
      end do
      
      do i=1, nvars
         if posvar(i) == 0 then
            vartype(i) = HYPRE_SSTRUCT_VARIABLE_CELL
         end if
         if posvar(i) == 1 then
            vartype(i) = HYPRE_SSTRUCT_VARIABLE_NODE
         end if

		 !Will replace this with new version of HYPRE, or avoid using face-variables
         if posvar(i) == 2 then
            vartype(i) = HYPRE_SSTRUCT_VARIABLE_XFACE
         end if
          if posvar(i) == 3 then
            vartype(i) = HYPRE_SSTRUCT_VARIABLE_YFACE
         end if
          if posvar(i) == 4 then
            vartype(i) = HYPRE_SSTRUCT_VARIABLE_ZFACE
         end if
      
      end do

      do i=1,ndim
         corner(1, i) = 1
      end do

      m = 1
	  
      do i=1,ndim
        do j=1,m
	       do k=1, ndim
              corner(j + m, k) = corner(m, k)
	          if k == i then 
			     corner(j + m, k) = corner(m, k) + blocksize(i) - 1
              end if	          
           end do
        end do
        m = 2 * m		
      end do
      
      do i = 1,2**ndim
         do j = 1,ndim !Looping through corners
            do k = ndim
               offcorner(i,j,k) = 0
               if corner(i, j) == blocksize(j) then 
			      offcorner(i,j,k) = blocksize(j) + 1
               end if
            end do
         end do
      end do
	  
	  !Lowest indexed corner on Face (i, low) = corner(1)
	  !Highest indexed corner on Face (i, low) = corner(2**ndim - 2**(i-1))
	  
	  !Lowest indexed corner on Face (i, high) = corner(1 + 2**(i-1))
	  !Highest indexed corner on Face (i, high) = corner(2**ndim)
	  
      do i=1, ndim
         connector(i, 1, 1, :) = corner(1, :)
         connector(i, 1, 2, :) = corner(2**ndim - 2**(i-1), :)
         connector(i, 2, 1, :) = corner(1 + 2**(i-1), :)
         connector(i, 2, 2, :) = corner(2**ndim, :)
		 
         off_connector(i, 1, 1, :) = offcorner(1, i, :)
         off_connector(i, 1, 2, :) = offcorner(2**ndim - 2**(i-1), i, :)
         off_connector(i, 2, 1, :) = offcorner(1 + 2**(i-1), i, :)
         off_connector(i, 2, 2, :) = offcorner(2**ndim, i, :)
      end do
	  
      nlref = (lrefine_min - lrefine_max) ** ndim !FIX THIS: WHERE ARE MAXIMUM AND MINIMUM REFINEMENT DEFINED???
      end subroutine hypre_init
	  
      subroutine terminate_hypre

      use sstruct_data
	  
      implicit none
	  
      call HYPRE_SStructFACDestroy2(solver)
      call HYPRE_SStructGridDestroy(hypre_grid);
      call HYPRE_SStructGraphDestroy(hypre_graph);
      call HYPRE_SStructMatrixDestroy(full_ioperator);
      call HYPRE_SStructVectorDestroy(old_grid);
      call HYPRE_SStructVectorDestroy(new_grid);
	  
	  
      end subroutine terminate_hypre	
	  
      subroutine evolve(imp_part, omatrix, dt)

      use , ONLY :: change_high, change_low
	  
      implicit none
      
      integer :: rerun
      real :: max_change_all
      real :: imp_part, dt
      real*8, pointer :: omatrix(:, :, :,:)
	  
      call hypre_setup()
	  
      rerun = 1
      while rerun == 1 do
	  
	  !"grid" is the grid as of the beginnning of the time-step
	  !old_grid is the copy of that grid in the format needed for HYPRE
	  !new_grid is the final result
	  
	  !Initialize trial-grid and newgrid (use oldgrid to produce grid)
	  
	     rerun = 0
	     call explicit(dt * (1 - imp_part), omatrix)
		 
         call hypre_loadvars(dt * imp_part, omatrix)
         call hypre_run_solver()
		 
	     call get_change(max_change_all)
		 
	     if max_change_all > change_high then
		    rerun = 1
			dt = dt * (change_low + change_high) / (2 * max_change_all)
	     end if
		 
	     if max_change_all < change_low then
		    dt = dt * (change_low + change_high) / (2 * max_change_all)
	     end if
		 
      end do
      
      contains


	  
      subroutine hypre_setup()

      use sstruct_data
      use grid_data
	  
      implicit none
      integer :: i,j,k,l,lb,blockID, , connected_block, int_err, part
      integer, dimension(ndim), dim_map, dim_dir, local_index
      integer*8 :: full_ioperator, new_grid, old_grid !information defining non-zero part of implicit operator matrix
      real, dimension(ndim) :: physical_location
      integer, dimension(ndim) :: part_index, neigh_index !for mapping
      integer, dimension(ndim) :: base_index, neigh_base
      integer :: neigh_part, var_index, ierr, adjust
      logical :: same_part
      logical, dimension(2**ndim) :: done_block
      logical, dimension(ndim) :: neighref
      integer, dimension(ndim - 1) :: surface_factor, surface_index, &
         surface_dir
      integer*8, pointer :: matindex(:,:)
		 
      call HYPRE_SStructGridCreate(MPI_COMM_WORLD, ndim, nlref, &  !nlred = "nparts" in HYPRE
	     &hypre_grid)
      
	  
      do lb=1, size(blocklist)
         blockID = blocklist(lb)
		 
         part = 0
	     do i=1, ndim
		    !for rectangular grid, may change for hex-grid 
		    !(use trapezoids with different orientations to construct hexes, do that in geometry-module)
            dim_dir(i) = 1
            dim_map(i) = i
			
			
            part = lref(blockID, i) * lrefine_max ** (i - 1) + part
			
		    call grid_position(blockID, 1, 0, &
			   physical_location, local_index, extents(1,:))
 		    call grid_position(blockID, product(blocksize), 0, &
			   physical_location, local_index, extents(2,:))
			   
            do j=1, ndim
		       extents(1,j) = extents(1,j) / lrefine(lb, j)
		       extents(2,j) = extents(2,j) / lrefine(lb, j)
            end do
            
         end do
		 
         call HYPRE_SStructGridSetExtents(hypre_grid, part, &
		    extents(1,:), extents(2,:))
		 
      end do
      do part = 0, nlref - 1
         call HYPRE_SStructGridSetVariables(hypre_grid, part, nvars,& 
	        vartypes, int_err)
	        
	 
 !        do i=1,ndim !not done in FLASH, mirroring FLASH, sort of, for grid/graph-setup
  !          do j=1,2
	!		   connected_block = neighbour(1, lb, i, j)
	!		   if connected_block > -1 then
	!		      call HYPRE_SStructGridSetNeighborPart(hypre_grid,blockID,&
	!				 offconnector(i, j, 1), offconnector(i, j, 2), connected_block, &
	!				 connector(i, j, 1), connector(i, j, 2), connect_map, connect_dir)
	!		   end if
	!	    end do
     !    end do
	 
      
      end do
	  
      call HYPRE_SStructGridAssemble(hypre_grid)
	  
      call HYPRE_SStructGraphCreate(MPI_COMM_WORLD, &
	     hypre_grid, &hypre_graph)
      
	  !Gather refinement level information
      do i=1, size(lref,1)
	     call MPI_AllReduce (lref(blockID), lref(blockID), 1, & 
            MPI_INT, MPI_MAX, MPI_COMM_WORLD, error)
      end do
	  
	  !Set up graph to accept variables as FEM with appropriate sparsity
      allocate(matindex(size(grid, 2) * nvar * size(offset_list, 3), 2))
      m = 0
	  
      do lb=1,size(blocklist)
         blockID = blocklist(lb)
	     part = 0
	     do k=1, ndim
		    part += lref(blockID, k) * &
			         lrefine_max ** (k - 1)
	     end do
         HYPRE_SStructGraphSetFEM(hypre_graph, part)		 
		 	  !Set up the graphs (non-zero parts of the hypre_grid, needed to handle sparsity)
	 
	 !REPACE WITH INDICES GIVEN IN OMATRIX
	     do i=1, size(grid, 2) !grid(var,i,lb)
	        !Find global index of i
	        call grid_position(blockID, i, 0, &
               physical_location, local_index, global_index) 
	        do var = 1, nvar
	   	       rownum = var + global_index * nvar
		       do j = 1, size(offset_list, 3)
			      columnnum = rownum + offset_list(var, 1, j) * nvar
				  columnnum += offset_list(var, 2, j) - var
				  m += 1
				  matindex(m, 1) = rownum
				  matindex(m, 2) = columnnum
		       end do
	        end do
         end do
	     
	     call HYPRE SStructGraphSetFEMSparsity(hypre_graph, part, &
		    size(matindex, 1), matindex)
      end do 
      deallocate(matindex)
	  	  !!!GATHER OFF-PROCESOR INFORMATION AFTER AMR-STEP INTO GRID_DATA!!!
	  !Set up couplings
      do i=1,ndim !Looping over dimensions
	     do j=1,2 !Looping over low/high surface perpendicular to dimension
!		    if neighleaf(blockID,i,j,0)==neighleaf(blockID,i,j,2**(ndim-1)) & !Block near "top-left" is same as block near "bottom-right"
!			   !Must still check for difference in resolution in the direction of the 
!			   then 
!	           cycle !Same level of refinement, same Part
!		    end if
	        do l=1,2**(ndim-1) !Looping over corners of surface
	           connected_block = neighleaf(blockID, i, j, l)  !neigh(j, i, lb)
			!CREATE NEW VARIABLE TO STORE REFINEMENT-LEVEL OF NEIGHBOURS (FULL VECTOR) OR BOUNDARY-ID WHEN DOING AMR
			!MUST have dimensions (maxblocks, ndim, 2 (high / low), 2**(ndim - 1) (LEAF blocks on that face))
		       if connected_block < 0 then !Block is on boundary
			!higher half of numbers: lower resolution
			!lower half of numbers: higher resolution
			      
		          cycle
		       end if
		       part = 0
			   neigh_part = 0
               do k=1, ndim
			      part = part + lref(blockID, k) * &
			         lrefine_max ** (k - 1) 
				  neigh_part = neigh_part + lref(connected_block,k) * &
			         lrefine_max ** (k - 1)					 
               end
			
		       if part .eq. neigh_part then
			      cycle
		       end
			
			   !Get base-indices for coupling FINISH THIS PART, FIX GRID_POSITION
               call grid_position(blockID, 1, 0, &
			      physical_location, local_index, base_index)
		       !Get grid-position of other block from relation, refinement-level, and this one's physical position
		       do k=1, ndim
			      neigh_base(k) = base_index(k) + blocksize * (j-1) * (ndim-i+1)
				  neigh_base(k) = neigh_base(k) * 2**(neighref(k)-lref(blockID, k))
			   end
			   !Add to base-indices as appropriate, depending on difference in refinement
		       
		    !Loop through indices of each, set up coupling

		    m = 1
		    do k=1,ndim
			   if k /= i then
			      surface_dir(m) = k
				  m = m + 1
		       end if
		    end do
		    do k=1,ndim - 1
			   surface_factor(k) = 1
			   do m=k+1,ndim
			      if k /= i then
				     surface_factor(k) = surface_factor(k) * blocksize(m)
				  end if
			   end do
		    end do
			
 !For every value in this direction, do all combinations with not-yet-done directions
		    do m=1,ndim-1
		       surface_index(m) = 0
		    end do
		    do k=1,product(blocksize) / blocksize(i)
			   surface_index(ndim - 1) = surface_index(ndim - 1) + 1
			   do m=1,ndim-1
				  if surface_index(ndim - m) == blocksize(surface_dir(ndim-m)) then
				     surface_index(ndim - m) = 0
				     surface_index(ndim - m -1) = surface_index(ndim - m -1) + 1
				  end if
			   end do
			   adjust = 0
			   do m=1,ndim-1
			      adjust = adjust + surface_index(m) * surface_factor(m)
			   end do

			   part_index(k) = base_index(k) + adjust
		       if neighref(k) > lref(blockID, k) then !fine -> coarse boundary
				  neigh_index(k) = neigh_base(k) + adjust / 2
			   end
		       if lref(blockID, k) < neighref(k) then !coarse -> fine boundary
			      neigh_index(k) = neigh_base(k) + adjust * 2
			   end
			   if lref(blockID, k) == neighref(k) then !no change in this dimension
				  neigh_index(k) = neigh_base(k) + adjust
			   end
		    
			
               do var_index = 1, nvars
	              call HYPRE_SStructGraphAddEntries (hypre_graph, part, &
		             part_index,var_index, neigh_part,neigh_index, var_index)
               end do  
            end do  
         end do  
      end do  
      end do      
	  
      call HYPRE_SStructGraphAssemble(hypre_graph)
	  
	  !Create matrix-object for operator-matrix, in space defined by the graph
      call HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, hypre_graph, &
	     &full_ioperator)
		 
      call HYPRE_SStructMatrixSetObjectType(full_ioperator,HYPRE_PARCSR)
      call HYPRE_SStructMatrixInitialize(full_ioperator) !full implicit operator matrix ready to be filled
	  
	  !Create and set up the vector-objects for the old grid (end of previous time-step) and new grid (end of new time-step)
      call HYPRE_SStructVectorCreate(MPI_COMM_WORLD, hypre_grid, & 
	     &old_grid)
      call HYPRE_SStructVectorSetObjectType(old_grid, HYPRE_PARCSR)
      call HYPRE_SStructVectorInitialize(old_grid)
      
      call HYPRE_SStructVectorCreate(MPI_COMM_WORLD, hypre_grid, & 
	     &new_grid)
      call HYPRE_SStructVectorSetObjectType(new_grid, HYPRE_PARCSR)
      call HYPRE_SStructVectorInitialize(new_grid)
           
      end subroutine hypre_setup
	  
      subroutine hypre_loadvars(dt, omatrix,imp_part)
	  
      use sstruct_data
      use grid_data
      use mechanism_data, only :: offset_list
	  
      implicit none
	  
      !real, pointer, dimension(:,:) :: block_grid !Dimensions are variable-index, local position-index
      integer :: lb, i, j, k, var, g_index_1d, acting_index_1d
      integer, dimension(ndim) :: local_index, global_index
      real, dimension(ndim) :: physical_location
      real, dimension(size(offset_list, 3))  :: imatrix_elem
	  !Can make this allocatable for position-dependent offsets, for long-range effects / multigrid
	  !To do long-range / multi-grid, create another part and establish couplings accordingly
	  
	  !Set up variable-order to match HYPRE
	  
      do lb = 1, size(blocklist)
         blockID = blocklist(lb)
		 
	     part = 0
         do i=1, ndim
		    part = part + lref(blockID, i) * &
			   lrefine_max ** (i - 1)  
         end
		 !Modify this to reduce memory-requirement of omatrix
	     
         do i=1, product(blocksize)
		    call grid_position(blockID, i, data_type, &
			   physical_location, local_index, global_index)
            do k=1, nvars
			   m += 1 !CYCLES THROUGH VARIABLES AT EACH POINT AND THEN THROUGH POINTS
			   var = varorder(k)
			   HYPRE_SStructVectorAddFEMValues(old_grid, part, &
			      global_index, grid(var,i,lb))
			   HYPRE_SStructVectorAddFEMValues(new_grid, part, & !Set initial guess
			      global_index, grid(var,i,lb))
		    end do
	     end do
		 
	     m = 0
	     do i=1, size(omatrix, 2)
		    if (i == 1) then 
			   call grid_position(blockID, omatrix((blockID, i, 3), data_type, &
			      physical_location, local_index, global_index)				  

		    else if (omatrix(blockID,i,3) .ne. omatrix(blockID, i-1,3)) then 
			   call grid_position(blockID, omatrix((blockID, i, 3), data_type, &
			      physical_location, local_index, global_index)

		    end if
		    if (i == product(blocksize) .or. )  then !Check for end of column using sparsity information 
		       m = global_index * nvars  + omatrix(blockID, i, 1) !New Row number
		
		
		    else  !Adding on to current row
			   
			
		    end if
			!Can probably replace most of omatrix with a single sparsity matrix. It should be the same.
!		    n = (global_index + omatrix(blockID, i, 4)) * nvars &
!			        + omatrix(blockID, i, 2) !Column number not needed. Column number used is determined by sparsity pattern
			
		    do k=1, nvars
			   do j = 1, size(omatrix, 2)
			      if omatrix(blockID, i, 1) == k then
			   !REPLACE WITH omatrix(blockID, elemnum, 5)
			   !elemnum being number of nonzero points
			   !blockID being block number
			   !5 being {var1, var2, unwrapped affected location, offset number (sparsity pattern gien elsewhere), coefficient}
			   !Just cycle through serial element, get rid of offset-list
			      imatrix_elem(j) = -omatrix(blockID, var1, i,j) * imp_part * dt !EDIT THIS TO WRITE TO THE CORRECT ELEMENT!!!!
			      if offset_list(var, 1, j) == 0 .and. &
				     offset_list(var, 2, j) == 0 then
			         imatrix_elem(j) = imatrix_elem(j) + 2.0 !Assumes that the operator-matrix includes the identity-component
			      end if
				  HYPRE_SStructMatrixAddFEMValues(full_ioperator, part, &
				   m, imatrix_elem)
               end do
            end do
         end do

		   ! HYPRE SStructMatrixSetBoxValues(full_ioperator, part, low_index, &
		    !   high_index, var, !COMPLETE THIS ROUTINE-CALL FOR BETTER EFFICIENCY, FIND OUT HOW TO DO IT
      end do

	  !BOUNDARY-CONDITIONS FOR FINE-COARSE BOUNDARIES ALREADY HANDLED BY GRID-INDEX UNIFICATION AND COMMUNICATION (FILLING) IN GRID-SECTION
      
	  !Set boundary-values at boundaries of grid (periodic conditions already handled by AMR)
	  
      HYPRE_SStructVectorGather(new_grid) !Needed for default object-type for storing vector, and for including face/node data
      HYPRE_SStructVectorGather(old_grid)
	  
      end subroutine hypre_loadvars
	  
      subroutine hypre_run_solver()
      
      use sstruct_data
      use grid_data
      implicit none
	  
      integer :: iterations
      integer*8 :: solver !Default solver is GMRES, I think
      double precision array(*) :: norm
	  
      integer :: part, lb, blockID, vect_index
      integer, dimension(ndim) :: base_index, i_location, 
      real, dimension(ndim) :: p_location
      real, dimension(product(blocksize)) :: values
      real :: max_change_local, max_change_all
	  
	  
      call HYPRE_SStructFACSetup2(solver,full_matrix,new_grid,old_grid)
      call HYPRE_SStructFACSolve3(solver,full_matrix,new_grid,old_grid)
!Get performance data
      call HYPRE_SStructFACGetFinalRelativeResidualNorm(solver, norm)
	  
      call HYPRE SStructFACGetNumIterations (solver, iterations)
	  
!Get results back to grid
	  
      do lb=1, size(blocklist)
         blockID = blocklist(lb)
		 
		 !Get global index of block, part
	     call grid_position(blockID, 0, 0, &
            p_location, i_location, base_index)
		 
	     upper_index = base_index + product(blocksize)
		 
	     part = 0
	     vect_index = 0
	     do i=1, ndim !Get part and unwrapped vector-index FINISH THIS
		    part = part + lref(blockID, i) * &
			   lrefine_max ** (i - 1)
			   !By default in HYPRE, the vector varies in the lowest dimension fastest
	     end do
	 
	     do var=1, nvar
		    call HYPRE_SStructVectorGetBoxValues(new_grid, part,base_index,&
               upper_index, var, values(:), ierr)
	        do i=1, size(values)
               newgrid(BlockID, var, i) = values(i)
	        end do
	     end do
		 
      end do
     
      end subroutine hypre_run_solver

      subroutine get_change (max_change_all)
	  
      use grid_data
      use mpi_data
	  
      implicit none
	  
      integer, intent out :: max_change_all
      integer max_change_local, max_change_all, change
      integer i, var, lb, blockID
	  
      max_change_local = 0
	  
	  !Find the largest change on this processor
      do lb=1, size(blocklist)
         blockID = blocklist(lb)
	     do var=1, nvar
	        do i=1, product(blocksize)
		      change = abs((newgrid(BlockID, var, i) - &
			    oldgrid(BlockID, var, i))/oldgrid(BlockID, var, i))
		      if change > max_change_local then
			     max_change_local = change
		      end if
		    end do
	     end do
      end do
	  
	  !Find the largest change overall
	  
      call MPI_AllReduce (max_change_local, max_change_all, 1, & 
         MPI_REAL, MPI_MAX, MPI_COMM_WORLD, error)
      end subroutine get_change
	  
      end subroutine evolve

	  