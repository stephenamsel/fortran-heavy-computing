module mcmc

     pure real*8 function step(param,i, parammin, parammax, scalesig)
	!Produces two Gaussian-distirbuted random variables, from Numerical Recipes.
	!Can increase speed by switching to a subroutine and modifying two parameters per run
	implicit none      
	integer i, turn 
      	real*8 z, sigma, scalesig
      	real*8 parammin(*), parammax(*), param(*)
      	real r, q, w
 
     	!Generate Gaussian-distirbuted random numbers
      	sigma= (parammax(i) - parammin(i)) * scalesig
      	r = rand()
      	q = rand()
      	!w = r + sqrt( -2.d0 * log(r*r+q*q)/(r*r+q*q) ) !w is now a Gaussian based on r, q
      	r = sigma*sqrt(-2*log(r))*cos(2*3.14159*q)
      	q = tan(2*3.14159*q)*r !Finding some way to use q could reduce runtime
           if (param(i) + r > parammin(i) .and. param(i) + r < parammax(i)) then
           	step=param(i)+r
           else
           	step = param(i)
           end if

      	return
      end

#if 0 THEN
This is a Markov Chain Monte Carlo wrapper around a set of calculators, "calctophat" and "calcqc".
It is optimized for parameter-space exploration as described here: http://arxiv.org/pdf/hep-ph/0407039v1.pdf .
As noted on Page 5, the optimal acceptance-rate of new points in parameter-space is 23.38%, but it is typically more effective to seek a 1 in 4 acceptance-rate because that allows much faster adjustment of step-sizes.
It achieves this goal by modifying the step-size dynamically, raising it when it is deemed unlikely that the acceptance-rate should be so low if the acceptance-probability is at or above 25% and lowering it when it is deemed unlikely that the acceptance-probability is at or below 25%. Smaller steps, in any system with a continuous posterior distribution, correspond to probabilities of acceptance closer to 1. As a result, if the initial value of the step is close to the correct order of magnitude, a sufficient burn-in should produce an optimized Markov chain.

This example uses a quality criterion calculated by calcqc and combines it with another crierion, a tophat distribution. It targets a specific quantity for the calculator, qc = targ, exploring the space targ +/- targerr.

The wrapper calculates only the next point in the chain. It it meant to be run in a loop with other elements.
For example, one might want to output the value of param in the loop and run other anaylses on the parameters found before feeding in the next point.

MC calls step, calctophat, and calcqc
#endif       
      subroutine MC(param, parammin, parammax, oldqc, targ, targerr, scalesig, sigmin, thinning)
      	implicit none
      
	real*8, intent(in) :: parammin(*), parammax(*), targ, targerr
	integer, optional, intent(in) :: thinning
	real*8, optional, intent(in) :: sigmin
	real*8, optional, intent(inout) :: scalesig 
!Used to translate between the Cartesian parameter-space explored and the spherical random walk of MCMC
	real*8, intent(inout) :: param(*), oldqc
	logical :: tophat
      	integer :: i, j, k
      	integer :: acc, rej
      	real*8 :: qc, z, !q is the quality-criterion, oldqc is the quality criterion from the previous point
      	real*8 :: Lnew, Lold, prob, test, step, newparam(*)
	
      	
	if( .not. present(scalesig) ) scalesig = 0.2
	if( .not. present(sigmin) ) sigmin = 0.01
	if( .not. present(thinning) ) thinning = 25 
!These values are chosen so that a random walk at the default step-scale could reasonably cross the explored space.
  
	allocate(newparam(size(params)))

      	acc = 0
      	rej = 0
      	do j=1,thinning !thinning to reduce correlation between output points in parameter-space.
!A higher value of thinning cuts correlations in the results while a lower value increases speed.

		do i=1,numparams
		         newparam(i)= step(param, i, parammin, parammax, scalesig)
	        end do

		test = rand()
!Gets likelihood of point, given data for qc (as defined, Z(gaussian) is strictly positive)
		Lold = 1.d0 - Erf(((oldqc-targ)*(oldqc-targ))**(1.d0/2.d0)/(0.0034d0))

      		call calcqc(newparam, qc)

      		Lnew = 1.d0 - Erf(((qc-targ) * (qc-targ))**(1.d0/2.d0)/(targerr))
      		prob = Lnew/Lold
!"tophat" is for a top-hat distribution (reject any point that is not "tophat"). 
!The "Quality criterion" is ["tophat" (1 or 0)] * [L(point)]

		call calctophat(tophat, newparam)
		if (test < prob .and. tophat) then
			do i=1, numparams
				param(i) = newparam(i)
			end do
	                oldqc = qc
			acc = acc + 1
	                if (acc > 2 .and. scalesig < 1.0d0) scalesig = scalesig * 1.25d0
				z = ((qc-targ)*(qc-targ))**(1.d0/2.d0)/(targerr)
				rej = 0
!If prob = 0.25, P(3 consecutive accepted points) = 1/64 < P(1 sigma)
           	else
			rej = rej + 1
			if (rej > 6 .and. scalesig > sigmin) scalesig = scalesig * 0.75d0
			acc = 0
                end if
!If prob = 0.25, P(7 consecutive rejected points) ~ 13% < P(1 sigma)
!These are not exact inverses. 1.25 * .75 is not 1. This is deliberate. It allows the step-size to approach the optimal value regardless of what it is.
      	end if
     
        return
      end

end module mcmc
