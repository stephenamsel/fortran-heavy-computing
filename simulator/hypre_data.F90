      Module hypre_data
      use grid_data, only :: nvar, ndim
	  
      implicit none
	  
      integer, dimension (nvar), save:: fieldID(nvar), fieldsSize(nvar)
      integer, allocatable, dimension (:), save :: elemFieldIDs, &
         nodeFieldIDs
      integer, save :: interleave, nodeNFields, elemNFields, nElems
	  
      end Module hypre_data