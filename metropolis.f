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

c **************************************************************
c
c This file contains the subroutines:  metropolis
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c
c **************************************************************

      subroutine metropolis(eol,currtem,acepta)

c===============================================================
c== SUBROUTINE FOR METROPOLIS UPDATE OF CONFIGURATIONS 	  ==
c==									 	  ==
c== CALLS: energy,addang,(rand),					  ==
c== dummy (function provided as argument)				  ==
c===============================================================


      include 'INCL.H'
      
c      declaration of variables
       real*8 e,t_col
	real KE,delta_i,delta_l,KE1,delta1,KE1_inter,KE2_inter,KE_ch1
	real minhits,theta,delta,difhits
	dimension vlvr1(mxrs),vlvr2(mxrs),c_vlvr(mxrs),c2_vlvr(mxrs)
	real chil1_vlvr(mxrs),chil2_vlvr(mxrs)
	integer desp1,desp2,Temp,cruza
	real time_i,eps,c1_eol,c2_eol,nval,betta
	real delta_ch1,delta_ch2,Einter,Einter1,Einter2

    	theta=15	
	minhits=20	
	eps=0.5	
        nval=0.8
        betta=0.8	
	
	delta_l=rand()	
	time_i = secnds(0.0)
c================================
c== Get Proposal configuration with soft perturbation==
c================================
c	collision type

	t_col  = rand()
    	
    		
    	
	if (t_col.LT.1.0) then	! collision on walls 

c================================
c== Decomposition ===
c================================
		jv = idvr(1+int(nvr*rand()))  ! select var.
		vrol = vlvr(jv)  ! save old
        	e = rand()
        	dv = axvr(jv)*e
        	vlvr(jv) = addang(vrol,dv)

		enw = energy()
	        delta = eol
	        delta_i = enw 

c	we evaluate if the molecule remain stable during a number hits 
		if(abs(enw-eol).LE.eps) then
			numhits=numhits+1
			difhits=minhits-numhits

c   	we apply perturbation 		
			if (difhits.lt.theta) then
				desp1=-49+int(mxrs*rand()*2)
				vlvr1=cshift(vlvr,desp1)
				desp2=-49+int(mxrs*rand()*2)
				vlvr2=cshift(vlvr,desp2)
				c_vlvr=vlvr	
				vlvr = vlvr1   
				c1_eol = energy()
				vlvr = vlvr2
				c2_eol = energy()

c	we evaluate the best structure 	
				if (c1_eol.le.enw.and.c1_eol.le.c2_eol) then
					vlvr=vlvr1
					delta_i = energy()
				else if (c2_eol.le.enw.and.c2_eol.le.c1_eol) then
					vlvr=vlvr2
					delta_i = energy()	
				else
					vlvr=c_vlvr
                                       delta_i=energy()								
				endif
												
			endif

		else
			numhits=1
		endif
	
c================================
c== check acceptance criteria ===
c================================
	        if (delta+KE.GE.delta_i) then
	          eol = enw
	          acepta = acepta + 1.0d0
		KE=((delta+KE)-delta_i)*delta_l        
	        else
	          vlvr(jv) = vrol
	        endif


	else		! inter-melecular collisions 
		write(*,*) ' inter molecular collision'
		
		c2_vlvr=vlvr
		jv = idvr(1+int(nvr*rand()))  ! select var.
		vrol = vlvr(jv)  ! save old
       		e = rand()
	      	dv = axvr(jv)*e
	     	vlvr(jv) = addang(vrol,dv)

		
c		we generate two structure applying crossover

		cruza =1+int(mxrs*rand())

		do i=1,mxrs
        		if (i.le.cruza) then 					
				chil1_vlvr(i)=c2_vlvr(i)
				chil2_vlvr(i)=vlvr(i)
			else
				chil1_vlvr(i)=vlvr(i)
				chil2_vlvr(i)=c2_vlvr(i)
			endif	          	
	       	enddo


c		synthesis criterion
		if (KE1_inter.ge.betta) then

			delta=energy()											
			c_vlvr=vlvr
			vlvr=c2_vlvr
			delta_i=energy()
			vlvr=chil1_vlvr	
			delta_ch1=energy()
			
			if (delta+delta_i+KE1_inter+KE2_inter.ge.delta_ch1) then
				KE_ch1= (delta+delta_i+KE1_inter+KE2_inter)-delta_ch1
				vlvr=chil1_vlvr
				delta=delta_ch1			
			else	

			endif				

		else

c		IntermolecularIneffectiveCollision

			delta=energy()											
			c_vlvr=vlvr
			vlvr=c2_vlvr
			delta_i=energy()
			vlvr=chil1_vlvr	
			delta_ch1=energy()
			vlvr=chil2_vlvr
			delta_ch2=energy()

			Einter1=delta+delta_i+KE1_inter+KE2_inter
			Einter2=delta_ch1+delta_ch2		
			Einter=Einter1-Einter2

				if (Einter1.ge.Einter2) then
					KE1_inter = Einter*delta_l
					KE2_inter = Einter*(1-delta_l)					

c					We evaluate the best structure					
					if (delta_ch1.le.delta_ch2) then
						vlvr=chil1_vlvr
						eol = energy()
														
					else
						vlvr=chil2_vlvr
						eol = energy()								
					endif
												
				else
					vlvr=c_vlvr
					eol=energy()	
				endif

		end if

	  endif			

      return
      end
