Program GPOPsim 
implicit none
	interface	
		subroutine popsim(mode,io_mode,nmk,nq,ngen,inh2,innq,selectlog)
		implicit none
			integer(1)			:: mode
			integer(1)			:: io_mode
			integer(4),optional :: nmk
			integer(4),optional :: nq
			integer(4),optional :: ngen
			integer(4),optional :: innq
			real(8),optional	:: inh2
			logical, optional	:: selectlog
		end subroutine
	end interface

	integer(1)			:: io_mode
	integer(1)			:: mode		 
	integer(4)			:: nmk
	integer(4)			:: nq
	integer(4)			:: ngen 
	integer(4)			:: io_log
	
	io_mode=1	! 1- interaction, 2- no interaction IN WINDOWS, 3/4 IN LINUX
	mode=1		! 1- main program, 2- subroutine
	io_log=39
	
	open(io_log, file='log.txt', status='replace')
	call popsim(mode,io_mode)
 	close(io_log)
end program


MODULE RANDOMM
IMPLICIT NONE
	INTEGER(4), SAVE	:: RSEED
	INTEGER(4)			:: seed1
	INTEGER(4)			:: seed2
	INTEGER(4)			:: seed3
CONTAINS

	FUNCTION randz( )
	IMPLICIT NONE       
		REAL(8)			:: randz

		RSEED=DMOD(16807D0 *RSEED, 2147483647D0)
		RANDZ=RSEED/2147483647.0
		RETURN
	END FUNCTION

SUBROUTINE startseed(x)
implicit none
	INTEGER(4)			:: rec(8)
	INTEGER(4)			:: x
	INTEGER(4)          :: clock    ! lixiujin
	
	if (x==0) then
	   call date_and_time (values=rec)
	   seed1=1.0E6*rec(8)+1.0E4*rec(3)+1.0E2*rec(6)+rec(7)
	  ! call system_clock(count=clock)
      ! seed1=clock !or seed1=clock+37      lixiujin added
	else 
		seed1=x
	end if
	rseed=seed1  !seed1-3 is random seed in pro. rseed is start seed.
	seed2=rseed**2
	seed3=rseed*(rseed-2)+5
END SUBROUTINE

Function rand1()
implicit none
	real(8)			:: rand1 	      
    real(8)			:: u1

       seed1=mod(seed1*171,30269)
       seed2=mod(seed2*172,30307)
       seed3=mod(seed3*170,30323)
       u1=seed1*1.d0/30269.d0+seed2*1d0/30307.d0+seed3*1.d0/30323.d0
       rand1=u1-int(u1)
END function

Function Randunif(a,b)
implicit none	
	integer(4)		:: a
	integer(4)		:: b
	integer(4)		:: randunif
	integer(4)		:: t
	real(8)			:: u

	u=randz()
	t=b-a+1
	Randunif=a+int(t*u)
	return
END function

Function Randunifreal(a,b)
implicit none	
	real(8)			:: a
	real(8)			:: b
	real(8)			:: randunifreal
	real(8)			:: t
	real(8)			:: u

	u=randz()
!	t=b-a+1         ! original code(error code)
    t=b-a           ! right code
	Randunifreal=a+t*u
	return
END function

	Function Stdnorm()
	implicit none
		real(8)			:: u1
		real(8)			:: u2
		real(8)			:: stdnorm
			
		U1=821.0 		  !801.0
		IF(U1.GE.254) GOTO 1100
		Stdnorm=U1
		U1=255
		RETURN
1100	u1=randz()
		U1=SQRT(-2*DLOG(U1))
		U2=6.2831852*randz( )
		stdnorm=U1*COS(U2)
		U1=U1*SIN(U2)
		RETURN
	END function

	Function gamma(a,b)
	implicit none
		real(4)			:: a
		real(4)			:: b
		real(8)			:: u1
		real(8)			:: u2
		real(8)			:: y
		real(8)			:: x
		real(8)			:: gamma
		real(8)			:: lamd
		real(8)			:: mu
		real(8)			:: cit
		real(8)			:: d
		real(8)			:: u
		real(8)			:: z
		real(8)			:: w
		integer(4)		:: i

		if (a>1.0E0 .and. amod(a,1.0E0)==0.0E0) then
			x=0.0E0
			do i=1,a
				x=log(randz())+x
			end do
			gamma=-1.0E0*x/b
			return
		elseif (a>1.0E0 .and. amod(a,1.0E0)/=0.0E0) then
			lamd=sqrt(2*a-1)
			mu=a**lamd
			cit=4.5
			d=1+log(cit)
	202		u1=randz()
			do i=1,50
				u2=randz()
			end do
			u=log(u1/(1.0E0-u1))/lamd
			x=a*exp(u)
			z=u1*u1*u2
			w=(a-lamd)*u+a-log(4.0)-x
			if (w+d-cit*z >=0.0E0 .or. log(z)<=w) then
				gamma=x/b
				return
			else 
				goto 202
			end if
		elseif (0.0D0 <= a .and. a <= 1.0D0) then 
	201		u1=randz()
			do i=1,50
			   u2=randz()
			end do
			y=b*u1
			if (y<=1) then
				x=y**(1.0E0/a)
			   if(u2<= exp(-x)) then
				    gamma=x/b
				    return
			   else
				    goto 201
			   end if			
		   elseif (y>1.0e0) then
				x=-1.0E0*dlog((b-y)/a)
			   if(u2<= x**(a-1.0E0)) then
				  gamma=x/b
				  return
			   else
				   goto 201
			   end if
		    else 
				call show(2 , 'QTL parameter error?',1)
		   endif
			
		else 
			call show(2 ,'RANDOM GAMMA DISTRIBUTION parameter error!',1)
			stop
		end if
	end function

	Function chisq(a)
	implicit none
		real(4)		:: a
		real(4)		:: chisq

		chisq=gamma( a/2.0E0, 0.5E0)

	END Function

function rexp(kk)     ! Exponential distribution
implicit none
	real(8)			:: kk
	real(8)			:: rexp

	rexp=-1.0D0*log(randz())/kk

end function


! Cholesky decomposition (test)
  subroutine Chde(V,n,L)
    implicit none
   integer,intent(in) :: n
   real(4),intent(in) :: V(n,n)
   real(4),intent(out) :: L(n,n)
   integer :: i,j,k
   real(4) :: sum1
      
   if(V(1,1)<=0.0E0) then
      write(*,*) "Cholesky decomposition failed for the first element no more than 0 "
   end if
           
    L(1,1)=sqrt(V(1,1))  ! determine numbers in col one
    do i=2,n
	  L(i,1)=V(i,1)/L(1,1)
    end do
   
    sum1=0.
	do j=2,(n-1)    !determine numbers in col j
	  do k=1,(j-1)
	    sum1=sum1+L(j,k)**2
	  end do

     if(V(j,j)-sum1<=0.0E0) then
         write(*,"(A,i,A,i,A)") "Cholesky decomposition failed for the",j,"*",j,"element no more than 0 "
      end if

     L(j,j)=sqrt(V(j,j)-sum1)
	 sum1=0.
	  do i=(j+1),n
	    do k=1,(j-1)
		  sum1=sum1+L(i,k)*L(j,k)
		end do
	
	   L(i,j)=(V(i,j)-sum1)/L(j,j)  
	   sum1=0.
     end do
   end do
   
   sum1=0.
   do k=1,(n-1)                !determine numbers in col n                  
     sum1=sum1+L(n,k)**2   
   end do
   if(V(n,n)-sum1<=0.0E0) then
         write(*,"(A,i,A,i,A)") "Cholesky decomposition failed for the",n,"*",n,"element no more than 0 "
   end if
   L(n,n)=sqrt(V(n,n)-sum1)

   return
end subroutine Chde 

subroutine multinorm(V,n,Y)
  implicit none
    integer,intent(in)  :: n
	real(4),intent(in)   :: V(n,n)
  !	real*8,intent(in)   :: u(n)
	real(4),intent(out)  :: Y(n)

	real(4),allocatable  :: L(:,:)
	real(4),allocatable  :: z(:)

	integer             :: i,j
	real(4)              :: u1

	allocate(L(n,n),z(n))
    
	z(1)=Stdnorm()
	do i=1,(n-1)
	  do j=1,50
	    u1=Stdnorm()
	  end do
	  z(i+1)=u1
    end do 

   call Chde(V,n,L)
    
  do i=1,n
    y(i)=0.0D0
    do j=1,i
	  y(i)=y(i)+L(i,j)*z(j)
	end do
 end do

 deallocate(L,z)

 return
end subroutine multinorm

END MODULE

subroutine printtime
	INTEGER(4)		:: DATE_TIME(8)
	CHARACTER(12)	:: REAL_CLOCK(3)
	integer(4)		:: io_log

	io_log=39
	CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                REAL_CLOCK (3), DATE_TIME)
	print *,real_clock(1)(7:8),'/',real_clock(1)(5:6),&
				'/',real_clock(1)(1:4),' ', &
			real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)
	write (io_log, *) real_clock(1)(7:8),'/',real_clock(1)(5:6),&
				'/',real_clock(1)(1:4),' ', &
			real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)

end subroutine

subroutine show(iomode, a, t)
implicit none
	integer(1)			:: iomode
	character*(*)		:: a
	integer(1),optional	:: t	! 1 output time, 0 no time show	
	INTEGER(4)			:: DATE_TIME(8)
	CHARACTER(12)		:: REAL_CLOCK(3)
	integer(4)			:: io_log
	
	io_log=39
	CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                REAL_CLOCK (3), DATE_TIME)
	write(io_log,'(a)') trim(a)	
	if (present(t) .and. t==1) then 
		write (io_log, *) real_clock(1)(7:8),'/',real_clock(1)(5:6),&
			'/',real_clock(1)(1:4),' ', &
			real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)
	endif
	if (iomode /= 4) then ! linux no interaction
		print *, trim(a)
		if ( present(t) .and. t==1) then 
			print *, real_clock(1)(7:8),'/',real_clock(1)(5:6),&
				'/',real_clock(1)(1:4),' ', &
				real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)		
		endif
	endif
	return
end subroutine

MODULE P_SIM
implicit none
	integer(4)	:: sim_nref

	integer(1), allocatable	:: qtype3(:,:)	! qtl genotype	
	integer(4), allocatable	:: pednew3(:)	! working vector for pedigree
	real(8), allocatable	:: vq3(:)		! qtl effect
	integer(4), allocatable	:: qtlid3(:)		! position for selected true qtl
	real(8), allocatable	:: bpv3(:,:)		! breeding values, phenotype values
END MODULE

SUBROUTINE POPSIM(mode,io_mode,nmk,nq,ngen,inh2,innq,l_select)
use p_sim
use randomm
implicit none
		! interface for show
	interface 
		subroutine show(iomode, a, t)
			integer(1)			:: iomode
			character*(*)		:: a
			integer(1),optional	:: t		
		end subroutine
	end interface
		! subroutine parameters
	integer(1)				:: mode			! 1-for main program,2- for subroutine 
	integer(1)				:: io_mode		! in/out put mode 1: interaction, 2: no interation, no screen print
	integer(4), optional	:: nmk			! number of markers
	integer(4), optional	:: nq			! number of qtl
	integer(4), optional	:: ngen			! number of generation
	integer(4), optional	:: innq			! no. qtl of input value 
	real(8),optional		:: inh2			! heritability of input value
	logical(4),optional		:: l_select 	!?? selection
	character(20)			:: parafile='para_in.txt'
		!POP_para
	integer(4)				:: psire		! total no. sire
	integer(4)				:: pdam			! total no. dam
	integer(4)				:: psize		! population size
	integer(4)				:: gsize		! genome size gsize=2*psize
	integer(4)				:: pgen			! total no. generations
	integer(4)				:: mutarule		! mutation rule:1- mutation; 0- no mutation
	integer(4)				:: materule		!?? mate rule: ??
	integer(4)				:: selerule		! selection rule:
	integer(4)				:: stage		! no. stage
	integer(4)				:: subnum		!?? sub number 
		!MQ_para
	real(8)					:: mratem		! mutation rate of marker
	real(8)					:: mrateq		! mutation rate of qtl
	real(4)					:: lenchr		! total length of chromsome
	integer(4)				:: nmark		! number of markers per chromsome
	integer(4)				:: nchr			! number of  chromsomes
	integer(4)				:: markrule		! rule for marker position
	integer(4)				:: tm			! total no. markers
	integer(4)				:: tq			! total no. putative qtl
	integer(4)				:: tqtl			! total no. qlt
	integer(4),SAVE			:: qtlrule		! rule for qtl 
	integer(4)				:: nqtl			! qtl number in output
	real(4),SAVE			:: qtlp1		! qtl parameter 1
	real(4),SAVE			:: qtlp2		! qtl parameter 2
!	real(4)					:: h2			! heritability
    real(4),allocatable	    :: h2(:)	    ! heritability
	real(4)					:: va			! Va
		!pub_matrix
	integer(4), allocatable	:: ped(:,:)		! pedigree 
	integer(4), allocatable	:: pednew(:,:)	! working vector for pedigree
	integer(4), allocatable	:: breds(:)		! selected sires
	integer(4), allocatable	:: bredd(:)		! selected dams
	integer(1), allocatable :: mtype(:,:)	! marker genotype
	integer(1), allocatable :: mtypenew(:,:)! working vector for marker genotype
	integer(1), allocatable :: mrec(:,:)	!?? record for marker mutation
	integer(1), allocatable	:: qtype(:,:)	! qtl genotype
	integer(1), allocatable :: qtypenew(:,:)! working vector for qtl genotype
	integer(1), allocatable :: qrec(:,:)	! record for qtl mutation
	integer(4), allocatable	:: qtlid(:)		! position for selected true qtl
	real(8), allocatable	:: rrate(:,:)	!?? recombination rate : 1 for m-q,2 for q-m 
	real(8), allocatable	:: vq(:,:)		! qtl effect,  
	real(8), allocatable	:: bpv(:,:,:)		! breeding values, phenotype values
		! test para
	integer(4)				:: recom		! counter of recombination
	integer(4)				:: mutm			! mutation for marker
	integer(4)				:: mutq			! mutation for qtl
	character(8)			:: timer1
	character(8)			:: timer2
	integer(4)				:: i,j,tgen,k
	logical					:: selog, l_count	
	character(40)			:: tmpchar

 !== lixiujin added===
   integer                  :: ntrait
   integer                  :: tqtl1,tqtl2,tqtl3
   integer                  :: ind_var     ! 0 to be continue trait ,1 to be threshold trait
   integer                  :: uncommon_QTL ! for 2 traits, 0 to be all QTL affecting two traits; 1 to be partial QTL
   real*4,allocatable       :: ve2(:,:)
   real*4,allocatable       :: ra(:),re(:)
   real*4                   :: qtleff_r1,qtleff_r2 ! r3
   real*4                   :: Tr   ! incidence of threshold trait
   real*4                   :: T    ! threshold value
   

	
!=============== Program start ======================
call show(io_mode , 'GPOPsim ___________Last modified 2014/10/17')

 
	if(.not. present(l_select)) i=1
	l_count=.false.
	if(l_count) open(99,file='hcount.out')

!=====STEP1=====parameter setting
	call time(timer1)			!timer1, record starting time
	call welcomeinf
	call readpara 

!=====STEP2=====base population installation
	call base_set

!=====STEP3=====MDE population 
	do i=1, pgen-10				
		call rand_mate(materule,psize)	
		call genotype()					
		call allele
		if(l_count)	call hcount(i)
		call select(selerule,psire,pdam)
		write(tmpchar,"(I5,3x,A20)") i,"generation finished!"
		call show(io_mode, tmpchar)
	end do

		mutarule=1	! no mutation
		call rand_mate(materule,psize)	
		call genotype()					
		if(l_count)	call hcount(i)
		call RECODEM(tm, mtype,gsize,mrec)
		call RECODEM(tq, qtype,gsize,qrec)
		if(l_count)	call hcount(i)
		call select(selerule,psire,pdam)
			write(tmpchar,"(I5,3x,A20)") i,"generation finished!"
			call show(io_mode, tmpchar)
	
	do i=pgen-9, pgen				
		call rand_mate(materule,psize)	
		call genotype()					
		if(l_count)	call hcount(i)
		call select(selerule,psire,pdam)
			write(tmpchar,"(I5,3x,A20)") i,"generation finished!"
			call show(io_mode, tmpchar)
	end do
	call show(io_mode,"Stage 1 finished!")

 

!=====STEP4=====one generation of pop - base population
	if (stage>1) then
		i=1
		call changepara(2) 
			call rand_mate(materule,psize)	
			call genotype()					
			if(l_count)	call hcount(i,1)
			if (mode==1) then
				call getqtl
				call outputmq
				call getBV
			   if(ind_var==1)  then        !  simulate threshold  
			     if(ntrait==1)   T=norminv(Tr,0.0,sqrt(1.0/h2(1)))
                 if(ntrait>=2)   T=norminv(Tr,0.0,sqrt(1.0/h2(1)))
					 write(*,*)"Threshold value=",T
			         call threshold
			    end if
				call show(io_mode, "Getbv passed.")
			endif		
			if (mode==2) then
				if(present(l_select))then
					if( l_select) then	!for gebv/phe selection, nqtl,h2 in para.in
					nq=tq
					call getgenotype(1)		!
						call show(io_mode, 'Getgenotype 1 passed')
						write(tmpchar, '(I12)') RSEED
						call show(io_mode, tmpchar) 
					call getqtl3 (qtlp1,qtlp2)		
						call show(io_mode, 'getqtl3 1 passed')
					call getbv3(1)	
					endif
				end if
			end if
			call select(selerule,psire,pdam)
			call output(pgen)
			i=1
			write(tmpchar,"(I5,3x,A20)") i,"generation finished!"
			call show(io_mode, tmpchar)
			call show(io_mode, "Stage 2 finished!")
	end if	
	tgen=1
!=====STEP5=====n generation of pedigree pop 
	if (stage>2) then
		do k=3, stage
		call changepara(k) 

			do i=1, pgen
				call rand_mate(materule,psize)	
				call genotype()					
				if(l_count)	call hcount(i+tgen)
				if (mode==1) then
					call getBV
                if(ind_var==1) call threshold  !lixiujin added to threshold trait
				endif
				if( mode==2) then
					if(present(l_select))then
						if(l_select) then
							call getgenotype(tgen+i)
							call getbv3(tgen+i)
						endif
					endif
				end if
				call select(selerule,psire,pdam)
				call output(tgen+i)
				write(tmpchar,"(I5,3x,A20)") tgen+i,"generation finished!"
				call show(io_mode, tmpchar)
			end do

		if (.false.) then
				do j=2, subnum
					seed1=rseed-1		!change random seed
					call baseback 		!setbase pop value
					tgen=tgen+(j-1)*pgen
					do i=1, pgen
						call rand_mate(materule,psize)	
							write(*,*) "rand_mate passed",tgen+i	
						call genotype()					
							write(*,*) "genotype passed ",tgen+i				
						if (mode==1) then
							call getBV
							write(*,*) "getbv passed    ",tgen+i
						endif
						call select(selerule,psire,pdam)
							write(*,*) "select passed   ",tgen+i
						call output(tgen+i)
							write(*,*) "putout passed   ",tgen+i
						write(*,"(I5,3x,A20)") tgen+i,"generation finished!"
					end do
					write (*,*) "Stage 3, subpop",j,"finished!"
				end do
		end if
		tgen=tgen+pgen
		write (tmpchar,'(A5, I4, A10)') "Stage",k," finished!"
		call show(io_mode, tmpchar)
		end do
	end if
	if (mode==2) then
		nmk=tm
		nq=tq
		ngen=tgen
	end if		
		
	call time(timer2)			!timer2, record endding time
	call overshow(timer1, timer2, seed1)
	if(l_count) close(99)

!======================================================================
CONTAINS
SUBROUTINE WELCOMEINF
	call show(io_mode, "          ========================================")
	call show(io_mode, "          |              GPOPSIM                 |")
	call show(io_mode, "          |Genomic POPulation SIMulation software|")
	call show(io_mode, "          |                 by Zhang Zhe 2010/03 |")
	call show(io_mode, "          |                  zhangzhecau@126.com |")
    call show(io_mode, "          |Function module-                      |") 
	call show(io_mode, "                       Multiple correlated traits|")             
    call show(io_mode, "          |                 by Xiujin Li 2014/10 |")
	call show(io_mode, "          |                lixiujin996@gmail.com |")
	call show(io_mode, "          ========================================")
END SUBROUTINE

! show runing information at the end of popsim
SUBROUTINE OVERSHOW(t1,t2,r)
implicit none
	character(*)		:: t1,t2
	integer(4)			:: r
	character(40)		:: tmpchar

	write(*,*)
	call show(io_mode, ' ')
	call show(io_mode,"        ========Congratulation========")
	call show(io_mode,"        Program end without any error!")
	write(tmpchar,"(9x,A8,"" --> "",A8)") t1,t2
	call show(io_mode, tmpchar)
	write(tmpchar, '(A26, I14)') "        With random seed  ", r
	call show(io_mode, tmpchar)
	call show(io_mode,"        ==============================")
END SUBROUTINE

SUBROUTINE READPARA 
implicit none
	integer(4)		:: x
	integer(4)		:: useed
     
    if(io_mode==1) then
         PRINT *, 'Please input the name of your parameter file.'
         read *, parafile
    else
	  call getarg(1,parafile)
      parafile=trim(parafile)
      
	end if
		open(11, file=parafile)
		    read (11,*) !line1
			read (11,*) ntrait,ind_var,Tr,uncommon_QTL
		
			if(ntrait==1) then
			    read (11,*) !line 2
				read (11,*) 
				read (11,*) !line 3
				read (11,*) 
            else
			   read (11,*) !line 2   
               allocate(ra(ntrait*(ntrait-1)/2))
			   read (11,*) (ra(j),j=1,ntrait*(ntrait-1)/2)
               read (11,*) !line 3
               allocate(re(ntrait*(ntrait-1)/2))
			   read (11,*) (re(j),j=1,ntrait*(ntrait-1)/2)
		    end if
			read (11,*) !line 4
			allocate(h2(ntrait))
			read (11,*) (h2(j),j=1,ntrait)
			read (11,*) !line 5
			read (11,*) nmark, nchr, markrule
			read (11,*) !line 6
			read (11,*) lenchr, mratem, mrateq 
			read (11,*) !line 7
			read (11,*) qtlrule, nqtl, qtlp1, qtlp2
			read (11,*) !line 8
			read (11,*) stage, subnum, useed 
			read (11,*) !line 9
			read (11,*) psire, pdam, psize, pgen 
			read (11,*) !line 10
			read (11,*) mutarule, materule, selerule 
		close (11)
	call startseed(useed)		!0-random seed, nonzero- fixed seed
	gsize=psize*2
	tm=nmark*nchr
	tq=(nmark-1)*nchr
	allocate(ped(3,psize))
	allocate(pednew(3,psize))
	allocate(breds(psire))
	allocate(bredd(pdam))
	allocate(mrec(tm,2))
	allocate(qrec(tq,2))
	allocate(mtype(tm,gsize))
	allocate(qtype(tq,gsize))
	allocate(mtypenew(tm,gsize))
	allocate(qtypenew(tq,gsize))
	allocate(rrate(2,tq))
	allocate(bpv(psize,2,ntrait))  ! lixiujin eidted
	mrec=0
	qrec=0 
	mrec(:,:)=2
	qrec(:,:)=2

END SUBROUTINE

SUBROUTINE changepara(X)
implicit none
	integer(4)	:: x
	integer(4)	:: i

	open (41, file=parafile)
	!do i=1,4+4*x
     do i=1,12+4*x  ! lixiujin added
		read (41,*)
	end do
		read (41,*) !line 7
		read (41,*) psire, pdam, psize, pgen 
		read (41,*) !line 8
		read (41,*) mutarule, materule, selerule 
	gsize=psize*2
	close(41)
END SUBROUTINE

SUBROUTINE READPARA3 
implicit none
	integer(4)		:: i,j,x
	integer(4)		:: useed

		open(11, file=parafile)
			read (11,*) !line 1
			read (11,*) nmark, nchr, markrule
			read (11,*) !line 2
			read (11,*) lenchr, mratem, mrateq 
			read (11,*) !line 3
			read (11,*) qtlrule, nqtl, qtlp1, qtlp2,h2 
			read (11,*) !line 4
			read (11,*) stage, subnum, useed 
			read (11,*) !line 5
			read (11,*) psire, pdam, psize, pgen 
			read (11,*) !line 6
			read (11,*) mutarule, materule, selerule 
			tgen=0
			do i=1, stage-1
				read (11,*)
				read (11,*) j, j, j, pgen 
				read (11,*) !line 6
				read (11,*) 
				tgen=pgen+tgen
			end do
		close (11)

	tq=(nmark-1)*nchr
END SUBROUTINE

subroutine base_set
use randomm
implicit none
	integer(1)	:: x1=1
	integer(1)	:: x2=2
	real(8)		:: dis
	real(8)		:: mean
	real(8)		:: kk1,kk2
	real(8)		:: pos
	integer(4)	:: i,j,loc,ct

!==step1== setup genotype of marker and qtl for all indivdiuals
	mtype=x1
	qtype=x1
!==step2== setup bred individual, pedigree
	deallocate(breds)
	deallocate(bredd)
	allocate(breds(psire))	
	allocate(bredd(pdam))
		!pedigree for the founder, no relationship
	do i=1,psize
		ped(:,i)=(/i,i,i/)		
	end do
	do i=1,psize/2
		breds(i)=2*i-1   ! odd number to be sire
		bredd(i)=2*i     ! even number to be dam
	end do
!==step3== setup position of all loci
	rrate=0.0D0
	select case (markrule)
	case (0)			!uniform distribution
		dis=0.5E0*lenchr/real(nmark-1)
		rrate=0.5E0*(1-exp((-2.0E0)*dis))			! Haldane's map function  dis=-ln(1-2r)/2
	case (1)			!test  exponent distri. for illmima 50k chip
		mean=lenchr/real(nmark-1)
!+++++++++++++++++++++++++++++++++++++++++++++++
		kk1=0.5*mean	!weight for mixed distribution
!+++++++++++++++++++++++++++++++++++++++++++++++
		kk2=mean-kk1
		ct=0
		open (98,file='mposition.out')
		do j=1,nchr  !get interval length
			do i=1,nmark-1
				loc=(j-1)*(nmark-1)+i
				rrate(1,loc)= Randunifreal(0.0d0, 2.0d0*kk1) + kk2*rexp(1.0D0)
			end do 
			dis=sum(rrate(1,loc-nmark+2:loc))	!total length of chrom j
			dis=lenchr/dis
			rrate(1,:)=rrate(1,:)*dis
			write (98,'(I6, I12, I12)') (j-1)*nmark+1, 0,0 !int(rrate(1,loc)*1.0D8) !bp
			do i=1, nmark-1						!rescale to chrom length and get map distance
				loc=(j-1)*(nmark-1)+i
				ct=ct+int(rrate(1,loc)*1.0D8)
				write (98,'(I6, I12, I12)') (j-1)*nmark+i+1, ct, int(rrate(1,loc)*1.0D8) !bp
				pos=randz()						!get qtl position in marker interval
				rrate(2,loc)=rrate(1,loc)*pos
				rrate(1,loc)=rrate(1,loc)-rrate(2,loc)
			end do
			ct=0
		end do
		rrate=0.5E0*(1-exp((-2.0E0)*rrate))				! Haldane's map function
		close(98)
	end select						
end subroutine

subroutine rand_mate(rule,size)
use randomm
implicit none
	integer(4), intent(in)		:: rule			! mating rule
	integer(4), intent(in)		:: size			! population size of next generation
	integer(4), allocatable		:: recs(:)		! times of mate for each sire
	integer(4), allocatable		:: recd(:)		! times of mate for each dam
	integer(4)					:: ts			! max no. dams to be mated to each sire
	integer(4)					:: td			! max no. offsprings to be produced for each dam
	integer(4)					:: maxs			! selected sires in last generation
	integer(4)					:: maxd			! selected dams in last generation
	integer(4)					:: id			! id number-- the id for new individuals
	integer(4)					:: xs			! random integer for sire
	integer(4)					:: xd			! random integer for dam
	integer(4)					:: step			! counter for error mating 
	logical						:: sib			! logical for sib or not
	logical						:: dif			! logical for same family or not
	integer(4)					:: i,j,k
	
	deallocate(pednew)
	allocate(pednew(3,size))
	maxs=ubound(breds,1)
	maxd=ubound(bredd,1)
	allocate(recd(maxd))
	allocate(recs(maxs))
1002	step=0
	recs=0
	recd=0

	dif=.true.
	ts=ceiling(real(maxd)/real(maxs))	 
	td=ceiling(real(size)/real(maxd)) 
	if (maxd< maxs) then
		call show(io_mode,'No. sire larger than dam. Program cannot run.',1)
		stop
	end if
	if (mod(real(size)/real(maxd),1.0E0) /=0.0E0)then
		call show(io_mode, 'litter size may not be equal for different families.')
		dif=.false.
	end if
	step=0	!timer for coded animal number.
	select case (rule)
	case (0)			!random mate
		do i=1,	maxs	 
			do j=1,ts
1021			xd=randunif(1,maxd)
				if (recd(xd) == 1) goto 1021			
				recd(xd)=1
				if (i==maxs .and. j==ts .and. dif==.false.) td=size-(maxd-1)*td
				do k=1, td
					step=step+1
					pednew(:,step)=(/step,breds(i),bredd(xd)/)		
				end do 
			end do		
		end do
	case (1)			!random but without sib mate  xxxxxxxxxxxxxxxxxxxxxxxxx 
		do i=1,size		 
1030	 	xs=randunif(1,maxs)
			if (recs(xs) >= ts) goto 1030
1031		xd=randunif(1,maxd)
			if (recd(xd) >= td) goto 1031
			sib=(ped(2,breds(xs))==ped(2,bredd(xd))) .or. (ped(3,breds(xs))==ped(3,bredd(xd)))
			if (sib) then !goto 1001		
!				write(*,*) "Sib mate"
				step=step+1
				if (step==100) then
					call show(io_mode,'Sib mate may not be avoided.',1)
					goto 1002
				end if
				goto 1031
			end if
			recs(xs)=recs(xs)+1
			recd(xd)=recd(xd)+1
			pednew(:,i)=(/i,breds(xs),bredd(xd)/)
		end do
	case (2)			! totally random mate
		do i=1,	size	 
 			pednew(1,i)=i
			pednew(2,i)=breds(randunif(1,maxs))
			pednew(3,i)=bredd(randunif(1,maxd))
		end do

	end select
	deallocate(ped)
	allocate(ped(3,size))
	ped=pednew
end subroutine			

subroutine select(rule,ns,nd)
use randomm
implicit none
	integer(4)				:: rule
	integer(4)				:: ns
	integer(4)				:: nd 
	integer(4),allocatable	:: rec(:)
	real(4), allocatable	:: tbv(:)
	integer(4)				:: size
	integer(4)				:: x
	integer(4)				:: tn		
	integer(4)				:: i,j

	deallocate (breds)
	deallocate (bredd)
	allocate(breds(ns))
	allocate(bredd(nd))
	breds=0
	bredd=0
	size=ubound(mtype,2)
	tn=size/4	!number of animals for one sex

	if (ns>tn .OR. nd>tn .OR. ns<=2 .OR. nd<=2) then   
		call show(io_mode, "Parameter error in select_module! code 1",1)
		stop
	end if

	allocate(rec(tn))
	rec=0
	select case(rule)
	case(0)			! random selection
		if (ns==tn ) then	!all sires be selected
			do i=1,tn
				breds(i)=2*i-1
!				bredd(i)=2*i
			end do
		else
			do i=1,ns
1100			x=randunif(1,tn)
				if (rec(x)==1) goto 1100		
				rec(x)=1
				breds(i)=2*x-1
			end do
			rec=0
		end if

		if (nd==tn ) then	!all dams be selected
			do i=1,tn
				bredd(i)=2*i
			end do
		else
			do i=1, nd 
1200			x=randunif(1,tn)
				if (rec(x)==1) goto 1200		
				rec(x)=1
				bredd(i)=2*x
			end do
		endif 

	case(1, 2)		!TBV/phe selection
		allocate(tbv(tn))
		if (nd==tn ) then	!all dams be selected
			do i=1,tn
				bredd(i)=2*i
			end do
		else				!should be select based on TBV
			if(rule==1)then
				do i=1,tn
					tbv(i)=bpv(2*i,1,1) !tbv
					rec(i)=2*i
				end do
			else
				do i=1,tn
					tbv(i)=bpv(2*i,2,1) !phe
					rec(i)=2*i
				end do
			endif
			call sort(tn, tbv, rec)
			do i=1,nd 
				bredd(i)=rec(tn+1-i)
			end do
		endif 

		if (ns==tn ) then	!all sires be selected
			do i=1,tn
				breds(i)=2*i-1
			end do
		else				!sires should be select based on TBV
			if(rule==1)then
				do i=1,tn
					tbv(i)=bpv(2*i-1,1,1) !tbv
					rec(i)=2*i-1
				end do
			else
				do i=1,tn
					tbv(i)=bpv(2*i-1,2,1) !phe
					rec(i)=2*i-1
				end do
			endif
			call sort(tn, tbv, rec)
			do i=1, ns 
				breds(i)=rec(tn+1-i)
			end do
		endif 
	case default
		call show(io_mode, "Parameter Error in select_module! code 2",1)
		stop
	end select
END subroutine

subroutine genotype
use randomm
implicit none
	integer(4)		:: pos(10)
	integer(4)		:: mloc
	integer(4)		:: qloc
	integer(4)		:: m
	integer(4)		:: q			!loc=location/position, size=pop size
	integer(4)		:: sid
	integer(4)		:: did
	integer(4)		:: oid
	integer(4)		:: id			!sire id/dam id
	integer(1)		:: ofip(2)
	integer(4)		:: i,j,k,newmut

	if (psize /= ubound(ped,2)) then
		call show(io_mode, "Parameter Error in genotype_module! code 1",1)
		stop
	end if

	deallocate(mtypenew) 
	deallocate(qtypenew) 
	allocate(mtypenew(tm,gsize))
	allocate(qtypenew(tq,gsize))
	
	if (mutarule==1) then			! NOT take mutation into consideration
		mratem=0.0D0
		mrateq=0.0D0
	end if
	do i=1,gsize				!animal i in new generation	
		oid=(i+1)/2
		sid=ped(2,oid)*2-2		!id in gamete matrix +1/2
		did=ped(3,oid)*2-2
		if (mod(i,2)==0) then 
			id=did
		else 
			id=sid
		end if
		do j=1, nchr
			m=(j-1)*nmark
			q=(j-1)*(nmark-1)
			ofip=(/randunif(1,2),0/)
			call loctype(ofip,0.5D0,mratem)			!first loci	
			if (ofip(2)==0 ) then
				mtypenew(m+1,i)=mtype(m+1,id+ofip(1))
			else  
				mtypenew(m+1,i)=mrec(m+1,1)
				newmut=mrec(m+1,1)
				do while (any(mtype(m+1,:) == newmut+1))
					newmut=newmut+1
				end do
				mrec(m+1,1)=newmut+1
			endif
			do k=2,nmark
				mloc=m+k
				qloc=q+k-1
!				call loctype(ofip,0.0E0,mrateq)				
				call loctype(ofip,rrate(1,qloc),mrateq)	
				if (ofip(2)==0 ) then
					qtypenew(qloc,i)=qtype(qloc,id+ofip(1))
				else  
					qtypenew(qloc,i)= qrec(qloc,1)
					qrec(qloc,1)=qrec(qloc,1)+1
					newmut=qrec(qloc,1)
					do while (any(qtype(qloc,:) == newmut+1))
						newmut=newmut+1
					end do
					qrec(qloc,1)=newmut+1
				endif
				call loctype(ofip,rrate(2,qloc),mratem)			
				if (ofip(2)==0 ) then  ! no mutation
					mtypenew(mloc,i)=mtype(mloc,id+ofip(1)) 
				else  ! mutation
					mtypenew(mloc,i)= mrec(mloc,1)
					newmut=mrec(mloc,1)
					do while (any(mtype(mloc,:) == newmut+1))
						newmut=newmut+1
					end do
					mrec(mloc,1)=newmut+1
				endif				
			end do
		end do
	end do
	deallocate(mtype)
	deallocate(qtype)
	allocate(mtype(tm,gsize))
	allocate(qtype(tq,gsize))
	mtype=mtypenew
	qtype=qtypenew
end subroutine

SUBROUTINE LOCTYPE(ip,rr,mur)
use randomm
implicit none
	integer(1)		:: ip(2)
	real(8)			:: rr
	real(8)			:: mur
	real(8)			:: u

	ip(2)=0
	u=randz()
	if (u<rr) then
		ip(1)=3-ip(1)
		recom=recom+1 !test
	end if
	u=randz()
	if (u< mur) then
		ip(2)=1
		mutm=mutm+1
	end if
	return
END SUBROUTINE

SUBROUTINE allele
implicit none
	integer(4)		:: maxa,mina,xt
	integer(4)		:: i,j,t

	do i=1,tm
		t=0;maxa=2;mina=2
		xt=max(mrec(i,1),mrec(i,2))
		do j=2,xt		
			if (any(mtype(i,:)==j)) then
				maxa=j
			elseif (t==0) then
				mina=j
				t=1
			end if			
		end do	 
		if (t==0) then	
			mrec(i,1)=maxa+1
			mrec(i,2)=maxa+1
		else
			mrec(i,1)=mina
			mrec(i,2)=maxa
		end if	
	end do
	
	do i=1,tq
		t=0;maxa=2;mina=2
		xt=max(qrec(i,1),qrec(i,2))
		do j=2,xt		
			if (any(qtype(i,:)==j)) then
				maxa=j
			elseif (t==0) then
				mina=j
				t=1
			end if			
		end do	 
		if (t==0) then	
			qrec(i,1)=maxa+1
			qrec(i,2)=maxa+1
		else
			qrec(i,1)=mina
			qrec(i,2)=maxa
		end if	
	end do
END SUBROUTINE

SUBROUTINE RECODEM(intm, inmtype,ingsize,inmrec)
use randomm
implicit none
	integer(4)		:: intm
	integer(4)		:: ingsize
	integer(1)		:: inmtype(:,:)
	integer(1)		:: inmrec(:,:)
	real(4)			:: dev
	real(4)			:: freq
	real(4)			:: minf
	integer(1)		:: recoderule
	integer(4)		:: i,j,c,lx
	
	recoderule=1
	do i=1, intm
	select case (recoderule)
	case (1) ! max freq. rule
!		if (mrec(i,2)<=2) cycle ci
			minf=1.0
!			do j=1, inmrec(i,2)		! 1 is not included in recoding process. NOT GOOD
			do j=2, inmrec(i,2)		! rule used before 20100622
				c=count(inmtype(i,:)==j)
				freq=real(c)/real(ingsize)
				dev=abs(0.5E0-freq)
				if (minf >= dev) then
					minf=dev
					lx=j
				endif
			end do
	case(2)	 ! accept all but min freq <0.05
		do
			j=inmrec(i,2)
			if (j < 3) then
				lx=1
				exit
			end if
			lx=randunif(1, j)
			c=count(inmtype(i,:)==lx)
			if (c==0) cycle
			freq=real(c)/real(ingsize)
			if (0.00001 < freq .and. freq <=0.05) then
				minf=randz()
				if (minf < 0.1) exit
			else 
				exit
			end if
	!		if (c/=0) exit
		end do
	end select
		do j=1,gsize
			if (inmtype(i,j)==lx) then
				inmtype(i,j)=1
			else
				inmtype(i,j)=2
			endif
		end do
	end do 
END SUBROUTINE

SUBROUTINE GETQTL
use randomm
implicit none
	real(4)					:: dev
	real(4)					:: freq
	real(4)					:: minf
	real(4),allocatable		:: tvg(:,:)
	real(4)             	:: stdvg
	integer(1)				:: qtl(tq)
	real(4)					:: u
	integer(4),allocatable	:: temp(:)
	character(40)			:: tmpchar
	integer                 :: cnqtl
	integer(4)				:: i,j,c,lx

	!==lixiujin added=====                for gamma distribution             
	 real*4,allocatable     :: vq1(:)
	 real*4,allocatable     :: V(:,:)
	 integer,allocatable    :: rec(:)
	 real*8                 :: x1,x2,x3,u1
	 real*4                 :: rq=0.0 ! the correlation of qtl effect
	 real*4,allocatable     :: corr_rq(:)
     real*4,allocatable     :: x(:,:)
	 integer                :: round,ncommon   
!	 real*4                 :: tvg_cov  
	 character(7),allocatable :: name(:)
	 real*4,allocatable       :: r_ra(:,:),r_re(:,:) 
	 integer                  :: k            
 
	
!	nqtl=int(qtlp1) !2009-5-16 nqtl read from parameter file dirctly
	qtl=0
	do i=1, tq  !recode qtl genotype
		if (any(qtype(i,:)/=qtype(i,1))) then
			minf=1.0E0   !recoding
			lx=1
			do j=1, qrec(i,2)
				c=count(qtype(i,:)==j)
				freq=real(c)/real(gsize)
				dev=abs(0.5E0-freq)
				if (minf >= dev) then
					minf=dev
					lx=j
				endif
			end do
		!	if (minf<=0.45E0) then		
				do j=1,gsize
					if (qtype(i,j)==lx) then
						qtype(i,j)=1
					else
						qtype(i,j)=2
					endif
				end do
				c=count(qtype(i,:)==j)
				freq=real(c)/real(gsize)
				qtl(i)=1
				tqtl=tqtl+1
!				if (abs(freq-0.5E0) <= 0.45E0) qtl(i)=1 !qtl position
		!	end if
		end if
	end do 
!	tqtl=sum(qtl)
	if (tqtl<=0) then
		CALL show(io_mode, "No QTL generated! Task failed! ")
		stop
	end if
	if (tqtl < nqtl) then
		call show(io_mode, 'tqtl litter than qtlp1')
		stop
	end if

	allocate(temp(tqtl))
	allocate(qtlid(nqtl))

	j=1
	do i=1,tq				!get selected qtl id in matrix
		if (qtl(i)==1) then
			temp(j)=i		! get qtl id
			j=j+1
		end if
	end do
    
	allocate(rec(tqtl))
	rec=0
	do i=1, nqtl
	100	lx=randunif(1,tqtl)
	  if(rec(lx)==1) goto 100
	  rec(lx)=rec(lx)+1
	  qtlid(i)=temp(lx)
	end do
    deallocate(rec)

	tqtl=nqtl
	allocate(vq(tqtl,ntrait))
    vq=0.0D0
   

   if(ntrait==2 .and.uncommon_QTL==1)then
	 if(ra(1)/=0.0E0) then

	   qtleff_r1=0.8;qtleff_r2=0.1
	 
	 else
	    
		 qtleff_r1=0.0
        ! qtleff_r2=Randunifreal(0,1)
		 qtleff_r2=0.5 
     end if
   end if

    round=0 
1003 select case (qtlrule)	!qtl effect      ! lixiujin added to simulate two trait
	case (0) !normal
	if(ntrait==1) then
	    stdvg=sqrt(qtlp2)
        do i=1,tqtl
			u=stdnorm()
			vq(i,1)=u*stdvg !+qtlp1
!			u=randunifreal(0.0E0,0.5E0)
!			if (u <0.5D0) vq(i)=-1*vq(i)
		end do
    endif
   
   if(ntrait==2 .and.uncommon_QTL==1)then 
     if(qtleff_r1/=0.0) then 
	   rq=ra(1)*(qtleff_r1+qtleff_r2)/qtleff_r1
	   if(abs(rq)>1.0) then
	      write(*,*) 'rq larger than 1.0 to adjust the setting value for genetic correlation!'
		  stop
	   end if
     end if
     tqtl1=0;tqtl2=0;tqtl3=0
     allocate(V(2,2),vq1(2))
	 V=0.0;vq1=0.0
     do i=1,tqtl
	  ! u1=randz()
      if(i<=qtleff_r1*tqtl) then
        V(1,1)=qtlp1
	    V(2,2)=qtlp2
	    V(1,2)=rq*sqrt(V(1,1)*V(2,2))
	    V(2,1)=V(1,2)
	    call multinorm(V,2,vq1)
        vq(i,1)=vq1(1)
		vq(i,2)=vq1(2)

		tqtl1=tqtl1+1
      elseif(i>qtleff_r1*tqtl.and.i<=(qtleff_r1+qtleff_r2)*tqtl)then
	     u=stdnorm()
	     vq(i,1)=u*sqrt(qtlp1)
         vq(i,2)=0.0D0
		 tqtl2=tqtl2+1
      else
	     vq(i,1)=0.0D0
	     u=stdnorm()
		 vq(i,2)=u*sqrt(qtlp2)
		 tqtl3=tqtl3+1
	  end if
    end do
   deallocate(V,vq1)

  end if
  
 if (ntrait>=2 .and.uncommon_QTL==0) then

  allocate(V(ntrait,ntrait),vq1(ntrait))
  V=0.0;vq1=0.0
  V(1,1)=qtlp1
  V(2,2)=qtlp2

  do i=3,ntrait
	 V(i,i)=qtlp2+real(i-2)
   end do
   k=0
   do i=1,ntrait
	  do j=(i+1),ntrait
	    k=k+1
	    V(i,j)=ra(k)*sqrt(V(i,i)*V(j,j)) 
		V(j,i)=V(i,j)
	  end do
    end do
    
	do i=1,tqtl
	    vq1=0.0
        call multinorm(V,ntrait,vq1)
		do j=1,ntrait
		  vq(i,j)=vq1(j)
		end do
	 end do
   deallocate(V,vq1)	   
  end if
  	     
!			u=randunifreal(0.0E0,0.5E0)
!			if (u <0.5D0) vq(i)=-1*vq(i)
  case (1) !gamma
    if(ntrait==1) then
       do i=1,tqtl
		vq(i,1)=gamma(qtlp1,qtlp2)
	    if(randz() <= 0.5E0) vq(i,1)=-1.0d0*vq(i,1)
	  end do
	end if

  if(ntrait==2 .and.uncommon_QTL==1)then 
     if(qtleff_r1/=0.0) then 
	   rq=ra(1)*(qtleff_r1+qtleff_r2)/qtleff_r1
	   if(abs(rq)>1.0) then
	      write(*,*) 'rq larger than 1.0 to adjust the setting value for genetic correlation!'
		  stop
	   end if
     end if

	 tqtl1=0;tqtl2=0;tqtl3=0
	 do i=1,tqtl
	   if(i<=qtleff_r1*tqtl) then
	 
	      x1=gamma(0.4*abs(rq),qtlp2)              ! x-Ga(a,b) E(x)=a/b var(x)=a/b**2
	      x2=gamma(0.4*(1.0-abs(rq)),qtlp2)
	      x3=gamma(0.4*(1.0-abs(rq)),qtlp2)
        vq(i,1)=x1+x2
		vq(i,2)=x1+x3
		u=randz()
        if(u>0.50E0) then
		  vq(i,1)=-1.0d0*vq(i,1)
  	      vq(i,2)=-1.0d0*vq(i,2)
		end if

	     if(rq<0) then
		    vq(i,1)=-1.0d0*vq(i,1)
		    vq(i,2)=vq(i,2)
            u=randz()
            if(u>0.50E0) then
		      vq(i,1)=vq(i,1)
		      vq(i,2)=-1.0d0*vq(i,2)
		    end if
	     endif
        
        tqtl1=tqtl1+1
       elseif(i>qtleff_r1*tqtl.and.i<=(qtleff_r1+qtleff_r2)*tqtl)then
	     vq(i,1)=gamma(0.4,qtlp2)
	     if(randz() <= 0.5E0) vq(i,1)=-1.0d0*vq(i,1)
         vq(i,2)=0.0D0
		 tqtl2=tqtl2+1
       else
	     vq(i,1)=0.0D0
		 vq(i,2)=gamma(0.4,qtlp2)
	     if(randz() <= 0.5E0) vq(i,2)=-1.0d0*vq(i,2)
		 tqtl3=tqtl3+1
	   end if
	 end do
     
   end if

   if (ntrait>=2.and.uncommon_QTL==0) then
      call show(io_mode,'please choose normal distribution because of No Mutiple-Gamma for 2 or above traits with uncommon_QTL==0!',1)
      write(*,*) "please choose normal distribution because of No Mutiple-Gamma for 2 or above traits with uncommon_QTL==0!"
      stop

   end if
		   
  case (2) !invers-chisquare
!		do i=1,tqtl
!			vq(i)=invchisq(??)
!		end do
 end select
! control rq 
  if(ntrait>=2 .and.(rq/=0.0 .or. uncommon_QTL==0)) then
    if(allocated(x)) deallocate(x)
	if(allocated(corr_rq)) deallocate(corr_rq)
    if(ntrait==2 .and.uncommon_QTL==1)then
	    ncommon=qtleff_r1*tqtl
    else 
	    ncommon=tqtl
	end if

    allocate(x(ncommon,ntrait),corr_rq(ntrait*(ntrait-1)/2))
	x=0.0;corr_rq=0.0
	do i=1,ntrait
	   x(:,i)=vq(1:ncommon,i)
	end do
    if(ntrait==2 .and.uncommon_QTL==1)then 
       call correl(x(:,1),x(:,2),ncommon,corr_rq(1))
       if(((rq-0.05)<=corr_rq(1)).and.(corr_rq(1)<=(rq+0.05))) then
             call show(io_mode,'producing the right correlation of qtl effect!',1)
	         write(*,*) rq,corr_rq(1),round
	         deallocate(x,corr_rq)
       else 
          round=round+1
	      if(round>=10000) then
	        call show(io_mode,'put new seed!',1)
		    stop
	      endif
          goto 1003
      end if   
   
  else 
   k=0
   do i=1,ntrait
	  do j=(i+1),ntrait
	    k=k+1 
	    call correl(x(:,i),x(:,j),ncommon,corr_rq(k))
        if(((ra(k)-0.05)<=corr_rq(k)).and.(corr_rq(k)<=(ra(k)+0.05))) then
            continue
        else 
            round=round+1
	        if(round>=10000) then
	           call show(io_mode,'put new seed!',1)
		       stop
	        endif
            goto 1003
        end if
	end do
  end do
  call show(io_mode,'producing the right correlation of qtl effect!',1)
  write(*,*) ra(1:(ntrait*(ntrait-1)/2)),corr_rq(1:(ntrait*(ntrait-1)/2)),round
  deallocate(x,corr_rq)
  	  
 end if	
end if
	
allocate(tvg(ntrait,ntrait))
tvg(:,:)=0.0E0
allocate(r_ra(ntrait,ntrait))
r_ra=1.0E0
allocate(ve2(ntrait,ntrait),r_re(ntrait,ntrait))
ve2=0.0E0;r_re=1.0E0


 if(ntrait==1) then
    do i=1,tqtl
	   freq=real(count(qtype(qtlid(i),:)==1))/real(gsize) !get qtl gene freq.
	   tvg(1,1)=tvg(1,1)+2.0E0*freq*(1.0E0-freq)*vq(i,1)*vq(i,1)
	end do
	vq(:,1)=sqrt(1.0E0/tvg(1,1))*vq(:,1)	!rescal qtl variance to 1
	tvg(1,1)=1.0E0
  end if
    ve2(1,1)=tvg(1,1)*(1.0E0/h2(1)-1.0E0)
 
   
 if(ntrait>=2) then
   do i=1,ntrait
	 do j=1,tqtl
		 freq=real(count(qtype(qtlid(j),:)==1))/real(gsize) !get qtl gene freq.
		 tvg(i,i)=tvg(i,i)+2.0E0*freq*(1.0E0-freq)*vq(j,i)*vq(j,i)		    
      end do
	     vq(:,i)=sqrt(real(i)/tvg(i,i))*vq(:,i)
		 tvg(i,i)=real(i)
   end do
   
   do i=1,ntrait
     do k=1,i-1
        do j=1, tqtl
	      freq=real(count(qtype(qtlid(j),:)==1))/real(gsize) !get qtl gene freq.
		  tvg(i,k)=tvg(i,k)+2.0E0*freq*(1.0E0-freq)*vq(j,i)*vq(j,k)
		 end do
		 tvg(k,i)=tvg(i,k)
		 r_ra(i,k)=tvg(i,k)/sqrt(tvg(i,i)*tvg(k,k)) 
		 r_ra(k,i)=r_ra(i,k) 
	 end do
   end do

     
	 do i=1,ntrait  
	    ve2(i,i)=tvg(i,i)*(1.0E0/h2(i)-1.0E0)
	 end do
	 k=0
     do i=1,ntrait
	    do j=(i+1),ntrait
		  k=k+1
		  ve2(i,j)=re(k)*sqrt(ve2(i,i)*ve2(j,j))
		  r_re(i,j)=re(k)
		  ve2(j,i)=ve2(i,j)
		  r_re(j,i)=r_re(i,j)
		  
	    end do
	 end do

   end if	  
  
     allocate(name(ntrait))
     do i=1,ntrait
	    write(name(i), "(A5,I2.2)") "trait",i
	 end do

	open(unit=21,file="Genetic_paramenters.out",status="replace")
       write(21,"('	',<ntrait>(A7,'	'))")(name(j),j=1,ntrait) 
	   write(21,"('Realized Additive Covariance matrix:')")
	   do i=1,ntrait
	      write(21,"(<ntrait>(f,'	'))") (tvg(i,j),j=1,ntrait)
	   end do
       write(21,"('Realized Additive Correlation matrix:')")
	   do i=1,ntrait
	      write(21,"(<ntrait>(f,'	'))") (r_ra(i,j),j=1,ntrait)
	   end do
       
       write(21,"('Residual Covariance matrix:')")
	   do i=1,ntrait
	      write(21,"(<ntrait>(f,'	'))") (ve2(i,j),j=1,ntrait)
	   end do
       write(21,"('Residual Correlation matrix:')")
	   do i=1,ntrait
	      write(21,"(<ntrait>(f,'	'))") (r_re(i,j),j=1,ntrait)
	   end do

	   write(21,"('h2:','	',<ntrait>(f,'	'))") (h2(j),j=1,ntrait)
	   
    close(21)
	 deallocate(name)
	 deallocate(r_re,r_ra)   	 	    
 
  	
	write (tmpchar, '(A19, I5)') "Total QTL number =",tqtl
	call show(io_mode, tmpchar)   
!	write (tmpchar, '(A19, 2x F12.8)') "Total QTL variance =", tvg(1)
	deallocate(tvg)

END SUBROUTINE


SUBROUTINE GETBV
use randomm
implicit none
	real(4)			    :: x
	real(4)             :: stdve
	real(4)			    :: u
	real(4),allocatable :: tbv(:)
	integer(4)		    :: genotype
	integer(4)	    	:: i,j,k

 !== lixiujin added===
	real*4,allocatable  :: e(:)
	

    if(allocated(tbv)) deallocate(tbv)
	allocate(tbv(ntrait))

    if(allocated(bpv)) deallocate(bpv)
	allocate(bpv(psize,2,ntrait))

	if(ntrait==1)  then
      stdve=sqrt(1.0E0/h2(1)-1.0E0)        ! residual variance =va(1.0/h2-1.0)
	  bpv=0.0E0
	  do i=1, psize
		tbv=0.0E0
		do j=1,tqtl		! tbv  phenotype
			genotype=qtype(qtlid(j),2*i-1)+qtype(qtlid(j),2*i)
			genotype =genotype - 3
			tbv(1)=tbv(1)+ vq(j,1)*genotype
!			select case (genotype)
!				case (2) !11
!					tbv=tbv+vq(j)
!				case (4) !22
!					tbv=tbv-vq(j)
!			end select
		end do
		bpv(i,1,1)=tbv(1)			! true breeding value
		u=stdnorm()
		bpv(i,2,1)=u*stdve+tbv(1)	! phenotype value
!		write(*,"(3F8.5)") tbv,u*stdve,bpv(i,2)
	 end do
   end if
  
   if(ntrait>=2) then                              !lixiujin added
	allocate(e(ntrait))
	e=0.0E0
	bpv=0.0E0

   	do i=1, psize
	  tbv=0.0E0
	  do k=1,ntrait
	   do j=1,tqtl		! tbv  phenotype
		   genotype=qtype(qtlid(j),2*i-1)+qtype(qtlid(j),2*i)
		   genotype =genotype - 3
		   tbv(k)=tbv(k)+ vq(j,k)*genotype
		end do
       end do
       
	   e=0.0E0
	   call multinorm(ve2,ntrait,e)
	   do k=1,ntrait
	     bpv(i,1,k)=tbv(k)			! true breeding value
		 bpv(i,2,k)=e(k)+tbv(k)
	   end do

     end do
	      
    !  bpv(i,2,1)=stdnorm()*sqrt(stdve(1))+tbv(1)
    !   do k=1,50
	!      u=stdnorm()                                 ! orign code to re=0
	!   end do
	!  bpv(i,2,2)=u*sqrt(stdve(2))+tbv(2)        
	 
	 deallocate(e)
  end if

END SUBROUTINE

Subroutine threshold            ! lixiujin added
   implicit none
    integer          :: i

     do i=1,psize
	     if( bpv(i,2,1)<T) then
		      bpv(i,2,1)=0.
		 else
		      bpv(i,2,1)=1. 
		 end if
	end do
end Subroutine threshold 
       

SUBROUTINE output(g)
!use randomm
implicit none
	integer(4)			:: bred
	integer(4)			:: ido
	integer(4)			:: idd
	integer(4)			:: ids
	character(12)		:: st
	integer(4)			:: i,j,g

	character(5),allocatable :: name(:)
	integer(4)          :: k

		write(st, "(I2.2, A9)") g, "MTYPE.out"
		open (20, file=st)
		write(st, "(I2.2, A9)") g, "QTYPE.out"
		open (21,file=st)
		write(st, "(I2.2, A6)") g, "BV.out"
		open (22,file=st)
		if(mode==1) then
		  if(ntrait==1) then
         	write (22, "(3A7,2A5, 2A10)") "ID", "Sire", "Dam", "Sex","Bred", "TBV", "Phe"
          end if

		  if(ntrait>=2) then
		    allocate(name(ntrait*2))
			do i=1,ntrait
			  write(name(2*i-1), "(A3,I2.2)") "TBV",i
              write(name(2*i), "(A3,I2.2)") "Phe",i
			end do
			    
    		write (22, "(3A7,2A5, <2*ntrait>(A10))") "ID", "Sire", "Dam", "Sex","Bred",(name(j),j=1,2*ntrait) 
			deallocate(name)
		  end if

		elseif(mode==2) then
			write (22, "(3A7,2A5)") "ID", "Sire", "Dam", "Sex","Bred"
		end if

		do i=1, gsize	!marker, qtl genotype output
			ido=int((i+1)/2)+g*10000
			write (20,"(2(2x I6.6),<tm>(I2))")  g,ido,mtype(:,i)   ! no, id, marker genotype
			if (mode==1) then
				write (21,"(2(2x I6.6),<tqtl>(I2))")  g,ido,(qtype(qtlid(j),i),j=1,tqtl)	! no, id, qtl genotype
			elseif(mode==2) then
				write (21,"(2(2x I6.6),<tq>(I2))")  g,ido, qtype(:,i)						! no, id, qtl genotype
			end if
		end do
		do i=1, psize	!breeding value, pedigree output
			bred=0
			if (mod(i,2)) then ! male
			    if (any(breds==i)) bred=1
			else				! female
			    if (any(bredd ==i)) bred=1
			end if
			if (g==0) then
				ido=ped(1,i)
				ids=0
				idd=0
			else 
				ido=g*10000+ped(1,i)
				ids=(g-1)*10000+ped(2,i)
				idd=(g-1)*10000+ped(3,i)
			end if
			if (mode==1) then
			  if(ntrait==1) then
                 	write (22,"(3(2x I6.6),2(I3),2(2x F18.12))") ido,ids,idd, &
			      & 2-mod(i, 2), bred, bpv(i,1,1), bpv(i,2,1)
             end if     
			 if(ntrait>=2) then
				   write (22,"(3(2x I6.6),2(I3),<2*ntrait>(2x F18.12))") ido,ids,idd, &
		     	 & 2-mod(i, 2), bred, ((bpv(i,k,j),k=1,2),j=1,ntrait)
			 end if

			elseif(mode==2) then
				write (22,"(3(2x I6.6),2(I3))") ido,ids,idd, 2-mod(i, 2), bred
			end if
		end do	
			
		close (20); close (21);close(22)
END SUBROUTINE

SUBROUTINE OUTPUTMQ
implicit none
	integer(4)			:: i,j
	character(len=10)	:: name
    
	open (22,file="MMAP.out")
	open (23,file="QMAP.out")
	write (23,'(<3>(I6))') nchr, nmark, nqtl
	if(ntrait==2 .and.uncommon_QTL==1)then 
	    write(23,*)'ratio of QTL:',tqtl1,'	',tqtl2,' ',tqtl3
    end if
	do i=1,tqtl		! qtl-name, id, effect
!		write (name,"(A3,I1,A1,I4.4)") "chr", int(qtlid(i)/(nmark-1))+1, &
!			& "-", mod(qtlid(i), (nmark-1))
!		write (23,"(A10,I5,2x F9.6)")  name, i, vq(i)!
      if(ntrait==1) then 
		write (23,'(I6, F18.12)') qtlid(i), vq(i,1)
	  end if
      if(ntrait>=2) then 
		write (23,'(I6,<ntrait>(F18.12))') qtlid(i), (vq(i,j),j=1,ntrait)
	  end if
   end do

	do i=1,nchr		! marker name, id
		do j=1,nmark
			write (name,"(A3,I1,A1,I4.4)") "chr",i,"-",j
			write (22,"(A10,I5)")  name, (i-1)*nmark+j	! 
		end do
	end do
END SUBROUTINE

subroutine Mloci()
implicit none
	logical::m(tm),q(tq)
	integer::hmark,hqtl,g
	integer::i,j
	
	m=0
	q=0
	hmark=0
	hqtl=0
	do i=1, tm
		m(i)=any(mtype(i,:)/=mtype(i,1))
	end do
	do i=1,tq
		q(i)=any(qtype(i,:)/=qtype(i,1))
	end do
	hmark=count(m)
	hqtl=count(q)
	write(*,"(A10,I5,1x,A19)") "There are ",hmark,"polymorphic markers."
	write(*,"(A10,I5,1x,A16)") "There are ",hqtl,"polymorphic qtls."
!	write(41,"(A5,A7,A5)") "Gen.","H-m","H-q"
!	write(41,"(I5,I7,I5)") g,hmark,hqtl

end subroutine

subroutine Hcount(gen,dd)
implicit none
	integer(4)		:: gen	!generation
	integer(4),optional	:: dd	!out put freq of current generation
	integer(4)		:: m1
	integer(4)		:: m2
	integer(4)		:: size 
	integer(4)		:: i,j,k
	real(4)			:: j1,j2,h,totalh,p,q


if(present(dd)) then
	open(42,file="HM.out")
!	open(43,file="Hq.dat")
	write(42,"(A6,2(A8))") "Pos", "Heter","Freq"  
!	write(43,"(A4,3(A8))") "ID", "Heter","F.1"  
endif
	size=ubound(qtype,2)
	totalh=0.0e0

	do i=1, tm
		h=0.0d0
		m1=count(mtype(i,:)==1)
		j1=real(m1)/real(size)
		do j=1, psize
			k=2*J-1
			if (mtype(i,k)/=mtype(i,k+1) ) h=h +1.0
		enddo
		h = h/real(psize)
		totalh=totalh+h
 		if (present(dd)) write(42,"(I6, F8.4, F8.4)") i,h,j1
	end do
	totalh=totalh/tm
!	if(present(dd))	write(42,"(A5,2x, F6.4)") "MeanH",totalh
!	write(99,"(A11,2x F6.4)") "MeanH mark=",totalh
	write(99,"(I6.4,2x F6.4)") gen ,totalh
	totalh=0
if (.false.) then
	do i=1,tq
		m1=count(qtype(i,:)==1)
		j1=real(m1)/real(size)
		if (j1 < 1.0E0) then
			j2=1-j1
			h=1-j1**2-j2**2
			totalh=totalh+h
			if (h/=0) write(43,"(I4,3(2x, F6.4))") i,h,j1,j2
		end if
	end do
	totalh=totalh/tq
	write(43,"(A5,2x, F6.4)") "MeanH",totalh
	write(*,"(A11,2x F6.4)") "MeanH qtl=",totalh
endif
if(present(dd))	close(42)
!	close(43)
end subroutine

SUBROUTINE baseback
implicit none
	integer(4)		:: i
	integer(4)		:: j
	integer(4)		:: js
	integer(4)		:: jd
	integer(4)		:: bred
	integer(4)		:: ido
	integer(4)		:: idd
	integer(4)		:: ids
	integer(4)		:: id
	character(18)	:: f20
	character(18)	:: f21
	character(18)	:: f22 
	character(12)	:: st
	integer(4)		:: nsire
	integer(4)		:: ndam
	integer(4)		:: nsize
	integer(4)		:: ngen

		open (19, file="para.in")
		do i=1,13
			read (19,*)
		end do
		read (19,*) nsire, ndam, nsize, ngen 
		close(19)

		deallocate(ped)
		deallocate(mtype)
		deallocate(qtype)
		deallocate(bredd)
		deallocate(breds)
		allocate(ped(3,nsize))
		allocate(mtype(tm,2*nsize))
		allocate(qtype(tq,2*nsize))
		allocate(bredd(ndam))
		allocate(breds(nsire))
		ped=0
		mtype=1
		qtype=1
		ndam=0
		nsire=0

		write(st, "(I2.2, A9)") ngen, "MTYPE.out"
		open (20, file=st)
		write(st, "(I2.2, A9)") ngen, "QTYPE.out"
		open (21,file=st)
		write(st, "(I2.2, A6)") ngen, "BV.out"
		open (22,file=st)

		read(22,*)
		do i=1, 2*nsize			!marker, qtl genotype input
			read (20,*)  id,id,mtype(:,i)   ! no, id, marker genotype
	!		write(*,*) i
	!		read (21,*)  id,id,	(qtype(qtlid(j),i),j=1,tqtl)	! no, id, qtl genotype
		end do
		do i=1, 2*nsize			!marker, qtl genotype input
	!		read (20,*)  id,id,mtype(:,i)   ! no, id, marker genotype
			read (21,*)  id,id,	(qtype(qtlid(j),i),j=1,tqtl)	! no, id, qtl genotype
		end do

		js=1
		jd=1
		do i=1, nsize	!pedigree input, bred animal input		
			read (22,"(3(2x I6.6),2(I3),2(2x F9.6))") ido,ids,idd, id, bred
			ids=mod(ids,10000)
			idd=mod(idd,10000)
			ped(:,i)=(/i,ids,idd/)
			if (bred==1) then
				if (mod(i,2)) then ! male
					bredd(js)=i
					js=js+1
				else				! female
					breds(jd)=i
					jd=jd+1
				end if
			end if
		end do	
			
		close (20); close (21);close(22)
END SUBROUTINE

SUBROUTINE sort(n,ra,rb)
implicit none
    real(4)			:: ra(*)
	real(4)			:: rra
	integer(4)		:: rb(*)
	integer(4)		:: l
	integer(4)		:: n
	integer(4)		:: ir
	integer(4)		:: rrb
	integer(4)		:: i,j

	l=n/2+1
	ir=n
	do
		if(l>1) then
			l=l-1
			rra=ra(l)
			rrb=rb(l)
		else
			rra=ra(ir)
			rrb=rb(ir)
			ra(ir)=ra(1)
			rb(ir)=rb(1)
			ir=ir-1
			if(ir==1) then
				ra(1)=rra
				rb(1)=rrb
				return
			endif
		endif
		i=l
		j=l+l
		do while(j<=ir) 
			if(j.lt.ir) then
				if(ra(j)<ra(j+1)) j=j+1
				endif
			if(rra<ra(j)) then
				ra(i)=ra(j)
				rb(i)=rb(j)
				i=j
				j=j+j
			else
				j=ir+1
			endif
		end do
		ra(i)=rra
		rb(i)=rrb
	end do
END SUBROUTINE sort

subroutine getgenotype(ge)
use p_sim
implicit none
	integer(4)			:: i,j,k
	integer(4)			:: ttsize
	integer(4)			:: ge
	integer(4)			:: io 
	character(40)			:: bvfile 
	character(40)			:: qfile 
	
	write (bvfile,'(I2.2,A6)') ge,'BV.out'
	write (qfile,'(I2.2,A9)') ge,'QTYPE.out'
	OPEN (101,file=bvfile) 
	OPEN (102,file=qfile)
	psize=0
	do 
		read(101,*,iostat=io)
		if (io/=0) exit
		psize=psize+1
	end do
	rewind(101)
	psize=psize-1	!for the first line
	gsize=psize*2

	if(allocated(pednew3))deallocate(pednew3)
	if(allocated(qtype3))deallocate(qtype3)
	allocate(pednew3(psize))
	allocate(qtype3(tq,gsize))
	pednew3=0
	qtype3=0
	k=0
	do i=1,psize
		read(102,*,iostat=io) j,pednew3(i),qtype3(:,2*i-1)
		if (io/=0) exit
		read(102,*,iostat=io) j,pednew3(i),qtype3(:,2*i)
		if (io/=0) exit
		k=k+1
	end do
	print *, k ,'individuals read finished'
end subroutine getgenotype

subroutine getqtl3(p1,p2)
use p_sim
use randomm
implicit none
	real(4)				:: p1,p2 !qlt parameter
	real(8)				:: tvg
	real(8)				:: stdvg
	integer(4)			:: i
	integer(4)			:: id 
	real(8)				:: u
	real(8),allocatable	:: freq(:)
	character			:: qmap*40
	character			:: tmpchar*40

	allocate(qtlid3(nqtl))
	allocate(vq3(nqtl))
	allocate(freq(nqtl))
	vq3=0.0D0
	freq=0.0d0
	qtlid3=0
	stdvg=1.0D0
	tvg=0.0D0

	i=1		!get qtlid
	do 
		id=randunif(1, tq)
		freq(i)=real(count(qtype3(id,:)==1))/real(gsize)
		if (freq(i)>=0.1D0 .and. freq(i)<=0.9D0 .and. .not.any(qtlid3==id)) then	
			qtlid3(i)=id
			i=i+1
		endif
		if (i > nqtl) exit !then 
		!	call show(io_mode, 'GetQTL3 error, program exit')
		!	stop
		!endif
	end do
	do i=1,nqtl	!sort id and freq, by id
		do j=i, nqtl
			if(qtlid3(i)>qtlid3(j)) then
				id=qtlid3(i)
				qtlid3(i)=qtlid3(j)
				qtlid3(j)=id
				u=freq(i)
				freq(i)=freq(j)
				freq(j)=u
			end if
		end do
	end do
	select case (qtlrule)	!get qtl effect
	case (0)	!normal
		do i=1,nqtl
			u=stdnorm()
			vq3(i)=u*stdvg 
		end do
	case (1)	!gamma
		do i=1,nqtl
			vq3(i)=gamma(p1,p2)
			if(randz() <= 0.5E0) vq3(i)=-1.0E0*vq3(i)
		end do
	case (2)	!invers-chisquare
		!do i=1,tqtl
		!	vq(i)=chisq(??)
		!end do
	end select
	tvg=0.0D0
	do i=1, nqtl
		tvg=tvg+2.0E0*freq(i)*(1.0E0-freq(i))*vq3(i)*vq3(i)
	end do
	vq3=sqrt(1.0D0/tvg)*vq3	!rescal qtl variance to 1
	
	write(tmpchar, '(A20, 2x I6)') "Total QTL number = ",nqtl
	call show(io_mode, tmpchar)
	write(tmpchar, '(A22, 2x F12.8)') "Total QTL variance = ", tvg
	call show(io_mode, tmpchar)
		!OUTPUT QMAP
	write(qmap,'(A8)')'QMAP.out'
	open(111,file=qmap)
	write(111,'(A)')'	POS		EFFECT		FREQ'
	do i=1,nqtl
		write (111,'(I5,F18.12,F10.4)') QTLID3(I),VQ3(I),FREQ(I)
	end do
	CLOSE(111)
end subroutine

SUBROUTINE GETBV3(gen)
use p_sim
use randomm
implicit none
	integer(1), allocatable	:: genotype(:)
	integer(1), allocatable	:: qqtype(:,:)
	real(8)			:: stdve
	real(8)			:: tbv
	integer(4)		:: i,j
	integer(4)		:: gen
	character		:: bvfile*40
		
!	stdve=sqrt(1.0E0/h2-1.0E0)	!assume vg=1.0
    stdve=0.0E0
!	if (allocated(qqtype)) deallocate(qqtype)
	if (allocated(bpv3)) deallocate(bpv3)
	allocate(genotype(nqtl))
	allocate(qqtype(nqtl,gsize))
	allocate(bpv3(psize,2))
	bpv3=0.0E0

	write(bvfile,'(I2.2,A7)') gen,'BV2.out'
	open(112,file=bvfile)
	write(bvfile,'(I2.2,A5)') gen,'Q.out'
	open(113,file=bvfile)
	write(112,*)'	ID			TBV			PHE'

	do i=1, nqtl
		qqtype(i,:)=qtype3(qtlid3(i),:)
	end do
	do i=1, psize
		tbv=0.0E0
		genotype=qqtype(:,2*i-1)+qqtype(:,2*i)-3
		do j=1,nqtl		
			tbv=tbv+ vq3(j)*genotype(j)
		end do
		bpv3(i,1)=tbv			! true breeding value
		bpv3(i,2)=stdnorm()*stdve+tbv	! phenotype value
		write(112,"(I6.6, 2F18.12)") pednew3(i),bpv3(i,:)
		write(113,'(I6.4, I8.6, <nqtl>(I2))') gen, pednew3(i),qqtype(:,2*i-1)
		write(113,'(I6.4, I8.6, <nqtl>(I2))') gen, pednew3(i),qqtype(:,2*i)
	end do
	close(112);close(113)

END SUBROUTINE


!=============lixiujin added norminv subroutine to calculate threshold value(simulate threshold trait)===================================================
Function norminv(p, mu, sigma)
   implicit none
   real*4     :: p       ! probability
   real*4     :: mu      ! mean
   real*4     :: sigma   ! standard deviation
   real*4     :: norminv
 
     norminv=(-sqrt(2.0)*sigma)*erfcinv(2.0*p)+mu

End Function


 Function erfcinv(y)
   implicit none
	real*4,intent(in)    :: y
	real*4               :: erfcinv
	
    integer              :: x=0 
    real*4               :: a(6)=(/1.37060048277853E-02, -0.30514157123572, 1.52430406921683, -3.05730326797099, 2.7104108320361, -0.886226926452692/)
    real*4               :: b(5)=(/-5.31993152326407E-02, 0.631194675226722, -2.43279656031073, 4.17508199298248, -3.32017038822143/)
    real*4               :: c(6)=(/5.50475133993694E-03, 0.227968721711412, 1.69759245777087, 1.80293316878195, -3.0933546798435, -2.07759567640438/)
    real*4               :: d(4)=(/7.78469570904146E-03, 0.32246712907004, 2.445134137143, 3.75440866190742/)
    real*4               :: ylow = 0.0485 
	real*4               :: yhigh = 1.9515

	real*4               :: q,r,u

    

    If(y >= ylow .And. y <= yhigh) Then

        q = y - 1
        r = q * q
        erfcinv =(((((a(1)*r+a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0)
        
    End If
    If (y > 0. .And. y < ylow) Then

        q = Sqrt(-2.0 * Log(y / 2.0))
        erfcinv = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0)
    End If

    If (y > yhigh .And. y < 2.0) Then

        q = Sqrt(-2.0 * Log(1.0 - y / 2.0))
        erfcinv = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0)

    End If

    u = (erfc(erfcinv) - y) / (-2.0 / Sqrt(3.1415926) * Exp(-erfcinv**2))
    erfcinv = erfcinv - u / (1.0 + erfcinv * u)
    
	if(y==0.)      erfcinv ="¡Þ"
    if(y==2.)      erfcinv ="-¡Þ"
	if(y<0.)       erfcinv = "NaN"
	if(y>2.)       erfcinv = "NaN" 

End Function

Function erfc(x)
  implicit none
    real*4    :: x
	real*4    ::erfc

    erfc = erfcore(x,1)
End Function

Function erfcore(x,jint)
  implicit none
  real*4,intent(in)   :: x
  integer,intent(in)  :: jint 
  real*4              :: erfcore
  real*4              :: xbreak = 0.46875
  
  real*4              :: a(5)=(/3.16112374387057, 113.86415415105, 377.485237685302, 3209.37758913847,0.185777706184603/)
  real*4              :: b(4)=(/23.6012909523441, 244.024637934444, 1282.61652607737, 2844.23683343917/)
  real*4              :: c(9)=(/0.56418849698867, 8.88314979438838, 66.1191906371416, 298.6351381974, 881.952221241769, 1712.04761263407, 2051.07837782607, 1230.339354798, 2.15311535474404E-08/)
  real*4              :: d(8)=(/15.7449261107098, 117.693950891312, 537.18110186201, 1621.38957456669, 3290.79923573346, 4362.61909014325, 3439.36767414372, 1230.33935480375/)
  real*4              :: p(6)=(/0.305326634961232, 0.360344899949804, 0.125781726111229, 1.60837851487423E-02, 6.58749161529838E-04, 1.63153871373021E-02/)
  real*4              :: q(5)=(/2.56852019228982, 1.87295284992346, 0.527905102951428, 6.05183413124413E-02, 2.33520497626869E-03/)

  real*4              :: y,z,xnum,xden,del

  integer             :: i
   
    If (Abs(x) <= xbreak) Then

        y = Abs(x)
        z = y * y
        xnum = a(5) * z
        xden = z
        do i = 1, 3
            xnum = (xnum + a(i)) * z
            xden = (xden + b(i)) * z
        end do
        erfcore = x * (xnum + a(4))/(xden + b(4))

        If (jint /= 0)   erfcore =1-erfcore
        If (jint == 2)   erfcore =Exp(z)*erfcore
    End If

    If (Abs(x) > xbreak .And. Abs(x) <= 4.0) Then

        y = Abs(x)
        xnum = c(9) * y
        xden = y
        do i = 1,7
            xnum = (xnum + c(i)) * y
            xden = (xden + d(i)) * y
        end do
        erfcore =(xnum + c(8))/(xden + d(8))
        If (jint /= 2) Then
            z = floor(y * 16.0) / 16.0
            del = (y - z) * (y + z)
            erfcore = Exp(-z * z) * Exp(-del) * erfcore
        End If
    End If

    If (Abs(x) >4.0) Then
       
        y = Abs(x)
        z = 1.0 / (y * y)
        xnum = p(6) * z
        xden = z
        do i = 1, 4
            xnum = (xnum + p(i)) * z
            xden = (xden + q(i)) * z
        end do

        erfcore = z * (xnum + p(5)) / (xden + q(5))
        erfcore = (1.0 / Sqrt(3.1415926) - erfcore) / y
        If (jint /= 2) Then
            z = floor(y * 16.0) / 16.0
            del = (y - z) * (y + z)
            erfcore = Exp(-z * z) * Exp(-del) * erfcore
            erfcore = 0
        End If
    End If
   
	if(jint==0) then
        If (x > xbreak)   erfcore = (0.5 - erfcore) + 0.5
        If (x < -xbreak)  erfcore = (-0.5 + erfcore) - 0.5
    elseif(jint==1) then
        If (x < -xbreak)  erfcore = 2.0 - erfcore
    else
        If (x < -xbreak) Then
            z = floor(x * 16.0) / 16.0
            del = (x - z) * (x + z)
            y = Exp(z * z) * Exp(del)
            erfcore = (y + y) - erfcore
        End If
    end if

End Function  


! calculate correlation between two samples
       subroutine correl(x,y,n,corr_r)
       real*4,intent(in)::x(:),y(:)
       real*4 corr_r
       integer n
       
       real*4 sum_xy,sum_x,sum_xx,sum_y,sum_yy,tmp_a
       sum_xy=0.0;sum_x=0.0;sum_xx=0.0;sum_y=0.0;sum_yy=0.0
       do i=1,n
          sum_xy= sum_xy+x(i)*y(i)
          sum_x = sum_x+x(i)
          sum_xx= sum_xx+x(i)**2
          sum_y=sum_y+y(i)
          sum_yy=sum_yy+y(i)**2                 
       enddo
       
       tmp_a=sqrt(n*sum_xx-sum_x**2)*sqrt(n*sum_yy-sum_y**2)
       corr_r=(n*sum_xy-sum_x*sum_y)/tmp_a
       end subroutine




!====================================================================================================

END SUBROUTINE POPSIM
