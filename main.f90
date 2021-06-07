        implicit none
        real, parameter :: N=1E8 ! total number of bacteria in large intestine
        
        character(22) :: filename,abc,ab,filename1
        
        !real,parameter  ::l=250.0  !Total length of large intestine (in cm) 
        real,parameter  ::l=1.0 
        real,parameter::c=0.008 !drift of protein in cm^2 per hour  
        integer, parameter :: t=6*3600 !Total time in seconds
        real, parameter ::dt=1.0
        !real, parameter :: dx=0.04
        real, parameter :: dx=1.0
        integer,parameter :: sp=int(l/dx)
        integer, parameter ::tp=int(t/dt)
        real, parameter :: bt=(N)/float(sp)
                !!!!!!!Giving controlled diet 86g of protein, 81g of fat and
                !373g of carbohydrates !!!!!!!!!!!!
        real,parameter :: p0=1.0! in hecto gram
        real,parameter :: f0=0.45 ! in hecto gram 
        real,parameter :: c0=1.0 ! in hectogram
        real,parameter :: fi0=2.0
 
        real:: D_b=0.0 !Difusion coefficient of bacteri
      
        integer,parameter :: ts1=int(3600.0/1.0)
                
                real,parameter :: delta2= 0.2E-12*0.25 !Nutrient consumption for host
                real,parameter :: delta= 0.2E-15*0.25  !Nutrient consumption per bacteria per second 
                real,parameter :: delta1=1.0E-12 ! mass of a bacteria
        real,parameter :: death=0.0006*0.005*0.3*0.0
        real,parameter :: K_s = 0.1 ! concentration of nutrient when the growth rate is half of its maximum value in mili molar.
        real ::cfl
        integer,parameter :: tnob=4,tnon=4
        real,parameter :: k = 1.0E08, k1=1.0E07, k3=1.0E06,k9=5.0E06,kl=50000000
        real,parameter :: ae=1000
        real, dimension(tnob,tnon) ::ds,pb,cb,fb,fib,ph,ch,fh,fih,pt,ct,ft,fit
        real,dimension(tnob,tnon) :: growth_p,growth_c,growth_f,growth_fi
        real, dimension(tnob,tnob) :: bs
        real,dimension(tnob,sp):: b,b_new
        integer :: i, j,ss,ts,xp1,xm1,un,m,um,q
        real, dimension(sp) :: x, c_p, c_f, c_c,c_fi,rn,ln, c_p_new, c_f_new,c_c_new,c_fi_new
        real ::asso,p_con,c_con,f_con,fi_con,a,asso1,asso2
        real,dimension(tnob) ::alpha_p,alpha_c,alpha_f,alpha_fi,deathb,growth,growthf,k0
        real :: bsum,t1
        un=1
        um=2
         cfl=c*(dt/dx)
        open(101, file='bd.dat',action='read',status='old')
        !open(102, file='test_inter.dat',action='read',status='old')

        open(102, file='bs.dat',action='read',status='old')
        open(12,file='test.txt',status ='new')
       !   open(12,file='rho_50np_21600.txt',status ='old')     
       open(103, file='k_values.dat',action='read',status='old')
     
      do ss=1,sp
       do i=1,tnob
         !if (i .eq. 1) then
         !b(i,ss)=(bt/tnob)
         !else
         b(i,ss)=(bt/tnob)
         !end if
         
       end do
      end do


     do i = 1, tnob
        read(101,*) (ds(i,j), j=1, tnon)

        read(103,*) k0(i)
        !write(*,*) k0(i)
     end do
        close(101)
     do i = 1, tnob
        read(102,*) (bs(i,j), j=1, tnob)
!        write(*,*) (bs(i,j), j=1,tnob)
     end do
       
        close(102)
do m=1,200
 !     write(1,*)m 
      ! q=m
      !write(ab,'(I0)') q
      do ss=1,sp
        x(ss)=ss*dx
       

        C_P(ss)=p0    !Concentration of Cnutrient in mili molar
        C_F(ss)=f0    !Concentration of nutrient Cin mili molar
        C_C(ss)=c0    !Concentration of nutrient Cin mili molar
        C_fi(ss)=fi0

      end do
   do ts=1,tp
! write(*,*) ts 
      do ss=1,sp
        xp1=rn(ss)
        xm1=ln(ss)
        bsum=0       
    do i=1,tnob
          pb(i,1) = delta*(c_p(ss)/(k_s+c_p(ss)))*(1.0/(exp(-ds(i,1)+1.0)))
          cb(i,2) = delta*(c_c(ss)/(k_s+c_c(ss)))*(1.0/(exp(-ds(i,2)+1.0)))
          fb(i,3) = delta*(c_f(ss)/(k_s+c_f(ss)))*(1.0/(exp(-ds(i,3)+1.0)))
          fib(i,4) =delta*(c_fi(ss)/(k_s+c_fi(ss)))*(1.0/(exp(-ds(i,4)+1.0)))

          alpha_p(i)=pb(i,1)/delta1
          alpha_c(i)=cb(i,2)/delta1
          alpha_f(i)=fb(i,3)/delta1
          alpha_fi(i)=fib(i,4)/delta1
         
          ph(i,1) = delta2*(c_p(ss)/(k_s+c_p(ss)))*(1.0/(exp(-ds(i,1)+1.0)))
          ch(i,2) =delta2*(c_c(ss)/(k_s+c_c(ss)))*(1.0/(exp(-ds(i,2)+1.0)))
          fh(i,3) =delta2*(c_f(ss)/(k_s+c_f(ss)))*(1.0/(exp(-ds(i,3)+1.0)))
          fih(i,4) =delta2*(c_fi(ss)/(k_s+c_fi(ss)))*(1.0/(exp(-ds(i,4)+1.0)))

          pt(i,1)  =   pb(i,1)+ph(i,1)
          ct(i,2)  =   cb(i,2)+ch(i,2)
          ft(i,3)  =   fb(i,3)+fh(i,3)
          fit(i,4) = fib(i,4)+fih(i,4)

          asso=0

         
         do j=1,tnob
          if (bs(i,j) < 0.0) then
!          write(*,*)m,ts,ss,i,j, bs(i,j)
          asso=asso+(bs(i,j))*b(j,ss)
         !asso=asso-b(j,ss) 

          end if
          !write(12,*)m,ts,ss,i,j, bs(i,j),b(j,ss)
          ! asso2=asso-b(j,ss)
          end do
          asso1=1
     
          !do j=1,tnob
          !if (bs(i,j) > 0.0) then
!          write(*,*)m,ts,ss,i,j, bs(i,j)
           !asso1=asso1*(bs(i,j)/(bs(i,j)+(b(j,ss)/bsum)))
           !end if
          !write(12,*)m,ts,ss,i,j, bs(i,j),b(j,ss)
           
         ! end do


 
          !write(*,*) i,k0(i)
!          growth_p(i,1)=alpha_p(i)*(asso1)*dt*(b(i,ss))*(1+asso/k0(i))
!
!          growth_c(i,2)=alpha_c(i)*(asso1)*dt*(b(i,ss))*(1+asso/k0(i))
          
!          growth_f(i,3)=alpha_f(i)*(asso1)*dt*(b(i,ss))*(1+asso/k0(i))

!          growth_fi(i,4)=alpha_fi(i)*(asso1)*dt*(b(i,ss))*(1+asso/k0(i))

           growth_p(i,1)=alpha_p(i)*(asso1)*dt*(b(i,ss))

          growth_c(i,2)=alpha_c(i)*(asso1)*dt*(b(i,ss))

          growth_f(i,3)=alpha_f(i)*(asso1)*dt*(b(i,ss))
          growth_fi(i,4)=alpha_fi(i)*(asso1)*dt*(b(i,ss))

          deathb(i) = death*dt*(b(i,ss))
          growth(i) =growth_p(i,1)+growth_c(i,2)+growth_f(i,3)+growth_fi(i,4)+(delta/delta1)*&
              &b(i,ss)*((asso/k)+(asso/kl)-((asso**2)/k*kl))
        
         growthf(i)=growth(i)
         b_new(i,SS)=b(i,ss)+growthf(i)-deathb(i)

       end do
            p_con=0.0
            c_con=0.0
            f_con=0.0
            fi_con=0.0

       do i=1,tnob
            p_con=p_con+(b(i,ss))*pt(i,1)
            c_con=c_con+(b(i,ss))*ct(i,2)
            f_con=f_con+(b(i,ss))*ft(i,3)
            fi_con=fi_con+(b(i,ss))*fit(i,4)
       end do

            C_P_new(ss)=C_P(ss)-p_con
            C_c_new(ss)=C_c(ss)-c_con
            C_f_new(ss)=c_f(ss)-f_con
            C_fi_new(ss)=c_fi(ss)-fi_con

!write(2,*) m,ts,ss,b(1,ss),b(2,ss),b(3,ss),b(4,ss), C_P_new(ss),C_c_new(ss),C_f_new(ss),C_fi_new(ss)
!write(2,*) m,ts,ss,b(1,ss),b(2,ss),b(3,ss),b(4,ss)
     end do
   do ss=1,sp   

        C_P(ss)=c_P_new(ss)
        C_C(ss)=c_C_new(ss)
        C_F(ss)=c_F_new(ss)
        C_Fi(ss)=c_Fi_new(ss)
        do i=1,tnob
          !  if (i .eq. 2) then
           !    b(i,ss)=0
           ! else  
               b(i,ss)=b_new(i,ss)
            !end if
               !do i=1,tnob
                 bsum=bsum+b(i,ss)   ! total number of bacteria at location ss
             !end do

        end do
   end do

50      format (12E11.4) 
        a=mod(ts,ts1)
        if(a .eq. 0) then
            write(abc,'(I0)') ts
            write(ab,'(I0)') m
            filename="rho_"//trim(ab)//"np_"//trim(abc)//".txt"
            un=un+1
            um=um+1
            open(unit=un,file=filename)
                do ss=1,sp
                 !  write(un,50)x(ss),b(1,ss),b(2,ss),b(3,ss),b(4,ss),b(5,ss),b(6,ss),b(7,ss),b(8,ss),b(9,ss),&
                 ! &b(10,ss),b(11,ss),b(12,ss),b(13,ss),b(14,ss),b(15,ss),b(16,ss),C_P(ss),c_c(ss),c_f(ss),c_fi(ss)
                !write(un,50) x(ss),b(1,ss),b(2,ss),b(3,ss),b(4,ss),b(5,ss),b(6,ss),b(7,ss),C_P(ss),c_c(ss),c_f(ss),c_fi(ss)

!                 write(un,50)b(1,ss)/bsum,b(2,ss)/bsum,b(3,ss)/bsum,b(4,ss)/bsum,b(5,ss)/bsum,b(6,ss)/bsum,&
!                        &b(7,ss)/bsum,b(8,ss)/bsum,b(9,ss)/bsum,b(10,ss)/bsum,b(11,ss)/bsum,b(12,ss)/bsum,&
!                        &b(13,ss)/bsum,b(14,ss)/bsum,b(15,ss)/bsum,b(16,ss)/bsum
!                write(un,50)b(1,ss)/bsum,b(2,ss)/bsum,b(3,ss)/bsum,b(4,ss)/bsum,b(5,ss)/bsum
               !  write(un,50)b(1,ss)/bsum,b(2,ss)/bsum        
              ! write(un,50)b(1,ss)/bsum,b(2,ss)/bsum,b(3,ss)/bsum   
               write(un,50)b(1,ss)/bsum,b(2,ss)/bsum,b(3,ss)/bsum,b(4,ss)/bsum
              
!            write(un,50) b(1,ss)/bsum,b(2,ss)/bsum,b(3,ss)/bsum,b(4,ss)/bsum,b(5,ss)/bsum,b(6,ss)/bsum,&
!                    &b(7,ss)/bsum,b(8,ss)/bsum,b(9,ss)/bsum,b(10,ss)/bsum,b(11,ss)/bsum,b(12,ss)/bsum

         !       write(un,50)b(1,ss),b(2,ss),b(3,ss),b(4,ss),b(5,ss),b(6,ss),b(7,ss),b(8,ss),b(9,ss),&
         !         &b(10,ss),b(11,ss),b(12,ss)

         end do
                 close(unit=un)
          end if
    end do
!write(5,*) m 
write(*,*) m
end do      
end
