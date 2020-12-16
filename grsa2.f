c Golden Ratio Simulated annealing 2 (GRSA2)
c Copyright (C) 2020  Dr. Juan Paulo Sánchez Hernández, and Dr. Juan Frausto Solis
c Copyright (C) 2005 Frank Eisenmenger, U.H.E. Hansmann, Shura Hayryan, Chin-Ku Hu

c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or (at
c your option) any later version.
c 
c This program is distributed in the hope that it will be useful, but
c WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c General Public License for more details.
c 
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
c USA.


c***********************GRSA************************************
c
c -This file contains a hybridization of simulated annealing called
c  Golden Ratio Simulated annealing 2 (GRSA2)
c -The GRSA2 has four fundamental parts: 1) The cooling scheme, 
c  2) The overheating strategies, 3) The convergence criteria, and 
c  4) soft pertubation based on CRO algorithm
c 
c -The firts version was developed by Dr. Juan Paulo Sánchez Hernández, 
c Fanny G. Maldonado Nava, and Dr. Juan Frausto Solis
c  2019-2020
c***************************************************************


      subroutine  grsa2(g_energy)

C --------------------------------------------------------------
C PURPOSE: SIMULATED ANNEALING SEARCH OF LOWEST-POTENTIAL-ENERGY
C          CONFORMATIONS OF PROTEINS
C
C CALLS: addang,energy,metropolis,outvar,outpdb,rgyr,setvar,zimmer
C
C ---------------------------------------------------------------
C ---------------------------------------------------------------

C------------------------------------------
C---------------------GRSA2-----------------
C------------------------------------------

      include 'INCL.H'

      character(8) x2
      logical lrand
      parameter(lrand=.true.)
	logical flag
	logical flag1
	logical flag2
	logical flag3
	logical flag4

        dimension vlvrm(mxvr)
	dimension idvrm(mxvr)
	dimension axvrm(mxvr)

        dimension vlvr_aux(mxvr)
	dimension idvr_aux(mxvr)
	dimension axvr_aux(mxvr)

	real*8 e,bangl
c	real inicio,final,total
	real initT,namefile
c	integer seed1,p,THREADS,ID,MaxThreads
 	

        flag = .true.
 	flag1 = .true.
 	flag2 = .true.
 	flag3 = .true.
 	flag4 = .true.
	s_window = 90.0d0
	t_flag=1

	s_slope = 0.000000000001d0
	slope = 0
	C = 0
	D = 0
	cnt = 1
	cnt2 = 0

	initT = secnds(0.0)
       nresi=irsml2(1)-irsml1(1) + 1     
c ==================
c == random start ==
c ==================
C== This start it is not necessary because we use a structure previously built
C       if(lrand) then
C        	do i=1,nvr
C         	iv=idvr(i)  ! provides index of non-fixed variable
C        	 e = rand(int(initT))
C        	 dv=axvr(iv)*e
C	         vr=addang(pi,dv)
C          	vlvr(iv)=vr
C        	enddo

C  	  end if		
        eol = energy()
	ymin = eol
 	epsilon = eol

c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
c      do i=1,ntlml
c        call outpdb(i,11)
c      end do

c==========================
c==	Simulated Annealing  ==
c==========================
      write(*,*) 'GOLDEN RATIO SIMULATED ANNEALING 2'
	call cpu_time(start)

c ====================================================
c == Initial parameters for Simulated Annealing ==
c ====================================================
	
c ====================================================
c ========== Parameters of 1iln3 ==============
c ====================================================	          

c====== This parameters were tunning by J. Frausto-Solís, H. Sanvicente-Sánchez, and F. Imperial-Valenzuela, "ANDYMARK: An Analytical Method to Establish Dynamically the Length of the Markov Chain in Simulated Annealing for the Satisfiability Problem," in Simulated Evolution and Learning, Springer Berlin Heidelberg, 2006, pp. 269–276.
c============================================================================================
		
	
	currtem = 44778080443048300000000000000000000000.0d0 !Initial Temperature of 1iln3
	stop_temp = 0.0000000200604823484395 !Final Temperature of 1iln3
      

	Tp0 = currtem*0.618
	Tp1 = Tp0*0.618
	Tp2 = Tp1*0.618
	Tp3 = Tp2*0.618
	Tp4 = Tp3*0.618

        alpha = 0.70
	blmax = 360.0d0
	bbeta = 1.00075046918563 

	
c =================================
c == Start the Simulated Annealing == 
c =================================
      do while(currtem.GE.stop_temp)
        propon = 0.0d0
        acepta = 0.0d0
        ycurr = 0.0d0
c =====================================
c == Start Metropolis ==
c == Crescent Markov Chain ==
c =====================================
	numhits=1
        do while(propon.LE.NINT(blmax))
          propon = propon + 1.0d0
          ycurr = eol
          
          if ((currtem.LE.TempEA) .AND. (ban1.EQ.1) ) then
	    count = count+1 
	    t_flag = 0    	  
	    write(*,*)'REHEAT', currtem,ymin
	    currtem = currtem+currtem
          end if

          call metropolis(eol,currtem,acepta)
c =================================================
c == Store the local lowest-energy conformation  ==
c =================================================
          if (eol.LT.ycurr) then
            ycurr = eol
          end if
c =================================================
c == Store the global lowest-energy conformation ==
c =================================================
          if (eol.LT.ymin) then
            ymin = eol
            write(*,*) 'MINIMA:', currtem, ymin

c ==================================
c == backup of Minimum point   ==
c == idvrm, axvrm, vlvrm  ==
c ==================================
       	do i=1,nvr
		  idvrm(i)  = idvr(i)
		  iv        = idvr(i)
		  axvrm(iv) = axvr(iv)
		  vlvrm(iv) = vlvr(iv)
       	enddo
          end if
c==================================================
c   The convergence criteria for dynamic equilibrium
c==================================================

			if ((currtem.LT.s_window).AND. (cnt.LT.3)) then
				D = D + ymin
				C = C + cnt * ymin
				cnt = cnt + 1
			end if

			if (cnt.EQ.3) then
		        slope = (((12 * C)-(6*(cnt-1)*D))/(cnt**3 - cnt))	
			cnt = 1	

			 if (slope.LT.s_slope) then
c----------------------------Overheating strategy-------------------------------------------
				currtem = currtem*0.618
				bangl= bangl+1
		         endif		
			
			  if (bangl.EQ.2) then
				currtem = stop_temp
			  endif

			endif


        enddo
c ========================
c == End of Metropolis ==
c ========================

c ==============================
c == Lower the Temperature ==
c ==============================
        currtem = alpha * currtem
	write(*,*) currtem,ymin
c ==============================
c Cooling Scheme by golden ratio
c ==============================

	if ((currtem.LT.Tp0).AND. (flag)) then
		alpha = 0.75d0
		flag = .false.
	else if ((currtem.LT.Tp1).AND. (flag1)) then
		alpha = 0.80d0
		flag1 = .false.
	else if ((currtem.LT.Tp2).AND. (flag2)) then
		alpha = 0.85d0
		flag2 = .false.
	else if ((currtem.LT.Tp3).AND. (flag3)) then
		alpha = 0.90d0
		flag3 = .false.
	else if ((currtem.LT.Tp4).AND. (flag4)) then
		alpha = 0.95d0
		flag4 = .false.
c--------------------Overheating strategy-------------------------------------------------
		currtem = currtem+currtem
	endif

c ==============================
c == Crescent Markov Chain ==
c ==============================
        blmax = bbeta * blmax

      enddo
	g_energy = ymin
c ====================================
c == End of the Simulated Annealing ==
c ====================================
	  namefile=g_energy*10000
	  write(x2,'(I8)') int(namefile) 
C	  open(17, file='./RESULTS/deltas'//x2//'.txt',status='new')
C	  write(17,*) 'Delta min y max',deltamin,deltamax	
          open(16, file='./RESULTS/Estruc'//x2//'.pdb',status='new')	  
C	  close(17)
c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
       do i=1,ntlml
         call outpdb(i,16)
       end do

       close(16)

      end
