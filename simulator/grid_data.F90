      module grid_data
	  
	     integer, save :: ndim
	     integer, dimension(ndim), save :: blocksize, indices_per_element
         integer, save :: block_indices, proc_indices
	     integer, dimension(maxrefine, ndim), save :: global_factor
	     integer, dimension(2, ndim), save :: block_boundary
	  
      end module grid_Data