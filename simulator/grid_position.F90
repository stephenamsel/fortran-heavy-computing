
      subroutine grid_position(blockID, data_index, data_type, &
         physical_location, index_location, global_index)

		use tree

		implicit none
!Finds physical location of point in grid, given index (cell-center)
!data_type: 0 for cell, 1 for corner, direction + 1 for cell-face

        integer, intent(in) :: blockID, data_type
        integer :: data_index, global_index
        integer, dimension(ndim) :: index_location
        real, dimension(ndim) :: physical_location
	
        integer i, indices_per_element(ndim), base_global_index(ndim)
	
        indices_per_element(ndim) = 1 !Number of indices per cell in that direction
		
        do i=1,ndim-1
           if (data_type == 0) indices_per_element(i) = indices_per_element(ndim + 1 - i) * blocksize(ndim - i)
           if (data_type == 1) indices_per_element(i) = indices_per_element(ndim + 1 - i) * (blocksize(ndim - i) + 1)
           if (data_type > 1) then
              if (ndim + 1 - i == data_type - 1) then
                 indices_per_element(i) = indices_per_element(ndim + 1 - i) * (blocksize(ndim - i) + 1)
              else
                 indices_per_element(i) = indices_per_element(ndim + 1 - i) * blocksize(ndim - i)
              end if
           end if
        end do

        do i=1,ndim
	       base_global_index(i)= int(coord(i,blockID)/dx(i, lrefine(blockID,i) + 0.5) !Handle rounding error
		   !TURN THIS INTO AN INTEGER TO MAKE IT ONE-DIMENSIONAL FOR HYPRE
        end do
		
        index_location(1) = data_index / indices_per_element(1)
        do i=2,ndim
           index_location(ndim) = mod((data_index / indices_per_element(i)), indices_per_element(i - 1))   
        end do
		
        do i=1,ndim
           if (data_type == 0) physical_location(i) = &
		      dx(refine(blockID, i) * index_location(i) + &
			     coord(i, blockID) + block%blockID%dx(i) / 2.0
           if (data_type == 1) physical_location(i) = &
		       dx(refine(blockID, i) * index_location(i) + coord(i, blockID)
           if (data_type > 1) then
              if (i == data_type - 2) then
                 physical_location(i) = dx(refine(blockID, i) * &
				    index_location(i) + coord(i, blockID)
              else
			   physical_location(i) = dx(i, refine(blockID, i)) * &
				  index_location(i) + coord(i, blockID) + dx(refine(blockID, i)/2.0
              end if
           end if
        end do
      end