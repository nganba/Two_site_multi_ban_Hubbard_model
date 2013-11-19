! two_site_states.f90
!
! Calculates the eigenvalues and eigenstates of 3 band two site hamiltonian for a given filling
! as a function of U, JH and spin-orbit coupling
!
! Author: Oinam Nganba Meetei
! Date:   08/02/13
!
!********************************************************************

program two_site_states
  use basis_states
  use hamiltonian
  use RXS
  implicit none
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  real(kind=double),parameter::pi=acos(-1.d0)
  real(kind=double)::U,JH,lambda
  integer::Nel,nbasis,i,j
  integer,dimension(1000)::basis
  real(kind=double),dimension(:),allocatable::eval,jz
  complex(kind=double),dimension(:,:),allocatable::hmat,jstate,jzmat
  complex(kind=double),dimension(:),allocatable::state
  real(kind=double)::Jmoment,Smoment,Lmoment,g_ratio,th,lmd,t,omega,L2,L3,theta,phi,norm
  real(kind=double),dimension(1000)::ene 
  integer,dimension(1000)::degen
  integer::itr,k
  real(kind=double),dimension(6)::occupation,occ
  character(10)::opt
  integer::deg
  complex(kind=double),dimension(3,3)::L2edge,L3edge
  real(kind=double),dimension(2,2)::L2polbasis,L3polbasis
  real(kind=double),dimension(3)::site_moment
  real(kind=double),dimension(3,3)::rot_x,rot_z

  write(*,*) 'U,Lambda,Nel'
  read(*,*) U,Lambda,Nel
  open(unit=1,file='L2edge.dat')
  open(unit=2,file='L3edge.dat')
  !U=40.d0
  JH=0.2*U
  !lambda=0.0d0
  !Nel=8
  !lmd=lambda
  !th=0d0
  !t=0d0
  degen=0
  ene=0
  call generate_basis(Nel,nbasis,basis)
  write(*,*) 'Number of basis = ',nbasis
  allocate(hmat(nbasis,nbasis),eval(nbasis),state(nbasis))
  do j=0,100
     !lambda=0.01*j
     U=0.5*j
     JH=0.2*U
     call generate_hamiltonian(U,JH,lambda,Nel,nbasis,basis,hmat)
     !call generate_hamiltonian(U,JH,lambda,Nel,nbasis,basis,hmat,th,lmd,t) 
     !do i=1,nbasis
     !   do j=1,nbasis
     !      if(hmat(i,j).ne.conjg(hmat(j,i))) then
     !         write(*,*) hmat(i,j)-conjg(hmat(j,i))
     !      endif
     !   enddo
     !enddo
     call diagonalize(hmat,eval)
     itr=1
     ene(1)=eval(1)
     k=1
     do i=2,nbasis
        if(abs(ene(itr)-eval(i)).gt.1d-5) then
           itr=itr+1
           ene(itr)=eval(i)
           degen(itr-1)=k
           k=1
        else
           k=k+1
        endif
     enddo
     ene=ene-ene(1)
     write(10,'(20(F9.5,2x))') U,lambda,ene(1:10)
     write(20,'(F9.5,2x,F9.5,2x,20(I4,2x))') U,lambda,degen(1:10)
     state=hmat(:,1)
     !write(*,*) 'degen', degen(1)
     allocate(jstate(nbasis,degen(1)),jzmat(degen(1),degen(1)),jz(degen(1)))
     jstate=hmat(:,1:degen(1))
     if(degen(1).gt.1) then
        call jz_matrix(jstate,nbasis,basis,jzmat)
        call diagonalize(jzmat,jz)
        do i=1,degen(1)
           hmat(:,i)=0d0
           do k=1,degen(1)
              hmat(:,i)=hmat(:,i)+jzmat(k,i)*jstate(:,k)
           enddo
        enddo
     endif
     !jstate=hmat(:,1:degen(1))
     !call jz_matrix(jstate,nbasis,basis,jzmat)
     !do i=1,degen(1)
     !   write(*,'(20(F8.5,2x))') jzmat(i,:)
     !enddo
     !write(*,*) degen(1)
     !Jmoment=j_moment(jstate,nbasis,basis)
     !Lmoment=l_moment(state,nbasis,basis)
     !Smoment=s_moment(state,nbasis,basis)
     !Jmoment=(-1.d0+sqrt(1.d0+4.d0*Jmoment))/2.d0
     !Lmoment=(-1.d0+sqrt(1.d0+4.d0*Lmoment))/2.d0
     !Smoment=(-1.d0+sqrt(1.d0+4.d0*Smoment))/2.d0
     !g_ratio= gyro_ratio(jstate,nbasis,basis)
     !write(*,*) g_ratio
     !write(30,'(20(F9.5,2x))') U,lambda,Jmoment,Lmoment,Smoment!,g_ratio*sqrt(abs(Jmoment))
     !write(*,'(20(F9.5,2x))') U,lambda,Jmoment,Lmoment,Smoment!,g_ratio*sqrt(abs(Jmoment))
     !write(*,'(20(F9.5,2x))') eval(1:15)
     !write(30, '(20(F9.5,2x))') U,lambda,Jmoment, g_ratio,(2.d0*Smoment-Lmoment)/Jmoment
     state=hmat(:,1)
     write(*,*) 'U,charge fluc',U, charge_fluctuation(state,nbasis,basis)/(Nel/2.d0)
     write(40,'(20(F9.5,2x))') U,lambda, charge_fluctuation(state,nbasis,basis)/(Nel/2.d0)
     occupation=0d0; deg=0; opt='j'
     do i=1,nbasis
        if(abs(eval(i)-eval(1)).gt.1d-5) exit
        state=hmat(:,i)
        call single_particle_occupation(state,nbasis,basis,occ,opt)
        occupation=occupation+occ
        deg=deg+1
     enddo
     occupation=occupation/deg
     write(50,'(20(F9.5,2x))') U,lambda,occupation
     deallocate(jstate,jzmat,jz)

     !state=hmat(:,3)
     !site_moment=siteJ_moment(state,nbasis,basis)
     !write(*,*) U,Lambda
     !write(*,*) 'Site moment'
     !write(*,'(3(F7.4,2x))') site_moment
     !norm=site_moment(1)**2 + site_moment(2)**2 + site_moment(3)**2
     !norm=sqrt(norm)
     !site_moment=site_moment/norm
     !theta=acos(site_moment(3))
     !if(abs(site_moment(1)).lt.1d-6) then
     !   phi=0d0
     !else
     !   phi=atan(site_moment(2)/site_moment(1))
     !endif
     !rot_z=0d0
     !rot_z(1,1)=cos(phi); rot_z(1,2)=sin(phi); rot_z(2,1)=-sin(phi); rot_z(2,2)=cos(phi); rot_z(3,3)=1.d0
     !site_moment=matmul(rot_z,site_moment)
     !write(*,'(3(F7.4,2x))') site_moment

     !call RXS_scattering_amplitude(U,JH,Lambda,Nel,L2edge,L3edge)
     !do i=1,3
     !   write(*,'(6(F9.5,2x))') L2edge(i,:)
     !enddo
     !write(*,*)
     !do i=1,3
     !   write(*,'(6(F9.5,2x))') L3edge(i,:)
     !enddo

     !do theta=0.d0,pi/2.d0,pi/80.d0
     !   L2polbasis=0d0; L3polbasis=0d0
     !   L2polbasis(1,1)=abs(L2edge(1,1))**2
     !   L2polbasis(1,2)=abs(sin(theta)*L2edge(1,2) + cos(theta)*L2edge(1,3))**2
     !   L2polbasis(2,1)=abs(-sin(theta)*L2edge(2,1) + cos(theta)*L2edge(3,1))**2
     !   L2polbasis(2,2)=abs(-(sin(theta)**2)*L2edge(2,2) + (cos(theta)**2)*L2edge(3,3) + sin(theta)*cos(theta)*(L2edge(3,2)-L2edge(2,3)))**2
     !   L3polbasis(1,1)=abs(L3edge(1,1))**2
     !   L3polbasis(1,2)=abs(sin(theta)*L3edge(1,2) + cos(theta)*L3edge(1,3))**2
     !   L3polbasis(2,1)=abs(-sin(theta)*L3edge(2,1) + cos(theta)*L3edge(3,1))**2
     !   L3polbasis(2,2)=abs(-(sin(theta)**2)*L3edge(2,2) + (cos(theta)**2)*L3edge(3,3) + sin(theta)*cos(theta)*(L3edge(3,2)-L3edge(2,3)))**2
     !   write(1,'(20(F10.7,2x))') U,Lambda,theta,L2polbasis(1,1),L2polbasis(1,2),L2polbasis(2,1),L2polbasis(2,2)
     !   write(2,'(20(F10.7,2x))') U,Lambda,theta,L3polbasis(1,1),L3polbasis(1,2),L3polbasis(2,1),L3polbasis(2,2)
     !enddo
     !write(1,'(20(F10.7,2x))') U,Lambda,abs(L2edge(1,1))**2!(abs(L2edge(1,1))**2+abs(L2edge(2,2))**2+abs(L2edge(3,3))**2)/3
     !write(2,'(20(F10.7,2x))') U,Lambda,abs(L3edge(1,1))**2!(abs(L3edge(1,1))**2+abs(L3edge(2,2))**2+abs(L3edge(3,3))**2)/3
  enddo

end program two_site_states
    
