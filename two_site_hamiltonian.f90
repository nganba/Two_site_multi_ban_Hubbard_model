! two_site_hamiltonian.f90
! 
! Generates the two site Hamiltonian including the effects of U, J_H and spin orbit 
! in a 3 band model. It is based on hmatgen_fermion.f90
!
! It also has subroutines to diagonalize the Hamiltonian (copied from hamiltonian.f90)
!
! Author: Oinam Nganba Meetei
! Date:   08/02/13
!
!******************************************************************************************

module hamiltonian
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer,parameter:: nsite=12
  integer:: lmax,ibin
  integer,save::j_init=-1
  integer,dimension(nsite)::lbit,nbitsps,nto2,nmltplt
  complex(kind=double),dimension(nsite,nsite),save::Hso,jx,jy,jz,sx,sy,sz,lx,ly,lz
  complex(kind=double),dimension(nsite,nsite),save::jxx,jyy,jzz,sxx,syy,szz,lxx,lyy,lzz
  complex(kind=double),dimension(:,:),allocatable,save::j2,s2,l2,l_z,s_z,j_z,j_x,j_y
  public:: generate_hamiltonian,diagonalize,j_moment,s_moment,l_moment,gyro_ratio, &
       charge_fluctuation,single_particle_occupation,jz_matrix,siteJ_moment

  interface diagonalize
     module procedure diagonalize_complex
     module procedure diagonalize_real
  end interface diagonalize

contains
  subroutine generate_hamiltonian(U,JH,lambda,Nel,nbasis,basis,hmat,th,lmd,t)
    real(kind=double),intent(in):: U,JH,lambda
    integer,intent(in)::Nel,nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:,:),intent(out)::hmat
    real(kind=double),intent(in),optional::th,lmd,t
    integer,dimension(nsite)::ksite,ksitenew,ksitenew2
    integer::id,jdg,i,j,k,idg,ires,izero,nbnd,nbonds,inew,inew2
    real(kind=double)::phi,phitot,theta,lambda2,t_hop
    call initialize(nbasis)
    if(size(hmat,1).ne.nbasis) then
       write(*,*) 'Hamiltonian array needs to have the same size a basis states'
       stop
    endif
    theta=0d0; lambda2=lambda; t_hop=1.d0
    if(present(th)) theta=th
    if(present(lmd)) lambda2=lmd
    if(present(t)) t_hop=t
    hmat=0d0
    izero=0
    !Single particle spin-orbit term
    ! Note: Basis used is in the order {dyz_up, dyzdn, dzx_up, dzx_dn, dxy_up, dxy_dn}
    Hso=0d0
    Hso(1,3)=-(0.d0,1.d0);  Hso(1,6)=1.d0
    Hso(2,4)=(0.d0,1.d0);   Hso(2,5)=-1.d0
    Hso(3,1)=(0.d0,1.d0);   Hso(3,6)=-(0.d0,1.d0)
    Hso(4,2)=-(0.d0,1.d0);  Hso(4,5)=-(0.d0,1.d0)
    Hso(5,2)=-1.d0;         Hso(5,4)=(0.d0,1.d0)
    Hso(6,1)=1.d0;          Hso(6,3)=(0.d0,1.d0)
    do jdg=1,nbasis
       id = basis(jdg)
       !Extract bit information
       call magxtr(id,ksite,izero)
       !write(*,'(20(I4,2x))') jdg,ksite 
       !Spin-orbit
       do i=1,6
          do j=1,6
             call cdag_i_c_j(ksite,ksitenew,inew,i,j)
             if(inew.eq.1) then
                call magcon(ksitenew,ires,izero)
                call binsrch(ires,idg,basis,nbasis)
                call phase(ksite,i,j,phi)            	 
                hmat(idg,jdg)= hmat(idg,jdg) - phi*(lambda/2.d0)*Hso(i,j)
             endif
          enddo
       enddo
       do i=7,12
          do j=7,12
             call cdag_i_c_j(ksite,ksitenew,inew,i,j)
             if(inew.eq.1) then
                call magcon(ksitenew,ires,izero)
                call binsrch(ires,idg,basis,nbasis)
                call phase(ksite,i,j,phi)            	 
                hmat(idg,jdg)= hmat(idg,jdg) - phi*(lambda2/2.d0)*Hso(i-6,j-6)
             endif
          enddo
       enddo
       !Crystal field
       !Tetragonal distortion
       hmat(jdg,jdg)=hmat(jdg,jdg) - 0.1d0*(ksite(5)+ksite(6)+ksite(11)+ksite(12))
       !Hopping term
       do i=1,6
          j=i+6
          call cdag_i_c_j(ksite,ksitenew,inew,i,j)
          if(inew.eq.1) then
             call magcon(ksitenew,ires,izero)
             call binsrch(ires,idg,basis,nbasis)
             call phase(ksite,i,j,phi)            	 
             hmat(idg,jdg)= hmat(idg,jdg) - phi*exp(-(0d0,1.d0)*theta)*t_hop
             !hmat(jdg,idg)= hmat(jdg,idg) - phi*exp( (0d0,1.d0)*theta)*t_hop 
          endif
          call cdag_i_c_j(ksite,ksitenew,inew,j,i)
          if(inew.eq.1) then
             call magcon(ksitenew,ires,izero)
             call binsrch(ires,idg,basis,nbasis)
             call phase(ksite,j,i,phi)            	 
             hmat(idg,jdg)= hmat(idg,jdg) - phi*exp(-(0d0,1.d0)*theta)*t_hop
             !hmat(jdg,idg)= hmat(jdg,idg) - phi*exp( (0d0,1.d0)*theta)*t_hop 
          endif
       enddo
       !Interaction
       !Intra-orbital
       do i=1,12,2
          j=i+1
          hmat(jdg,jdg)= hmat(jdg,jdg) + U*ksite(i)*ksite(j)
       enddo
       !Inter-orbital (same spin)
       !site 1
       do i=1,4
          do j=i+2,6,2
             hmat(jdg,jdg)= hmat(jdg,jdg) + (U-3.d0*JH)*ksite(i)*ksite(j)
          enddo
       enddo
       !Site 2
       do i=7,10
          do j=i+2,12,2
             hmat(jdg,jdg)= hmat(jdg,jdg) + (U-3.d0*JH)*ksite(i)*ksite(j)
          enddo
       enddo
       !Inter-orbital (opposite spin)
       !site 1
       do i=1,4
          k=i+3
          if(mod(i,2).eq.0) k=i+1
          do j=k,6,2
             hmat(jdg,jdg)= hmat(jdg,jdg) + (U-2.d0*JH)*ksite(i)*ksite(j)
          enddo
       enddo
       !site 2
       do i=7,10
          k=i+3
          if(mod(i,2).eq.0) k=i+1
          do j=k,12,2
             hmat(jdg,jdg)= hmat(jdg,jdg) + (U-2.d0*JH)*ksite(i)*ksite(j)
          enddo
       enddo
       !Spin-flip term
       !site 1
       do i=1,4,2
          do j=i+2,6,2
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,i,i+1) !cdag_up*c_dn
             if(inew.eq.1) then
                call phase(ksite,i,i+1,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,j+1,j) !cdag_dn*c_up
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,j+1,j,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) - phitot*JH
                endif
             endif
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,i+1,i) !cdag_up*c_dn
             if(inew.eq.1) then
                call phase(ksite,i+1,i,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,j,j+1) !cdag_dn*c_up
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,j,j+1,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) - phitot*JH
                endif
             endif
          enddo
       enddo
       !site 2
       do i=7,10,2
          do j=i+2,12,2
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,i,i+1) !cdag_up*c_dn
             if(inew.eq.1) then
                call phase(ksite,i,i+1,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,j+1,j) !cdag_dn*c_up
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,j+1,j,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) - phitot*JH
                endif
             endif
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,i+1,i) !cdag_up*c_dn
             if(inew.eq.1) then
                call phase(ksite,i+1,i,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,j,j+1) !cdag_dn*c_up
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,j,j+1,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) - phitot*JH
                endif
             endif
          enddo
       enddo
       !Pair Hopping term
       !site 1
       do i=1,4,2
          do j=i+2,6,2
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,j,i) !cdag_up*c_up
             if(inew.eq.1) then
                call phase(ksite,j,i,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,j+1,i+1) !cdag_dn*c_dn
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,j+1,i+1,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) + phitot*JH
                endif
             endif
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,i,j) !cdag_up*c_dn
             if(inew.eq.1) then
                call phase(ksite,i,j,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,i+1,j+1) !cdag_dn*c_up
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,i+1,j+1,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) + phitot*JH
                endif
             endif
          enddo
       enddo
       !site 2
       do i=7,10,2
          do j=i+2,12,2
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,j,i) !cdag_up*c_up
             if(inew.eq.1) then
                call phase(ksite,j,i,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,j+1,i+1) !cdag_dn*c_dn
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,j+1,i+1,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) + phitot*JH
                endif
             endif
             phitot=1.d0
             call cdag_i_c_j(ksite,ksitenew,inew,i,j) !cdag_up*c_dn
             if(inew.eq.1) then
                call phase(ksite,i,j,phi)
                phitot=phitot*phi
                call cdag_i_c_j(ksitenew,ksitenew2,inew2,i+1,j+1) !cdag_dn*c_up
                if(inew2.eq.1) then
                   call magcon(ksitenew2,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksitenew,i+1,j+1,phi)
                   phitot=phitot*phi
                   hmat(idg,jdg)= hmat(idg,jdg) + phitot*JH
                endif
             endif
          enddo
       enddo
    enddo
  end subroutine generate_hamiltonian



  subroutine single_particle_occupation(state,nbasis,basis,occupation,opt)
    implicit none
    complex(kind=double),dimension(:),intent(in)::state
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    real(kind=double),dimension(:),intent(out)::occupation
    character(10),intent(in),optional::opt
    complex(kind=double),dimension(6,6)::SO_evec
    real(kind=double),dimension(6)::SO_eval
    character(10)::option
    integer,dimension(nsite)::ksite,ksitenew
    integer::i,j,idg,jdg,id,INFO,orb,ires,izero
    real(kind=double)::phi
    complex(kind=double),dimension(6)::occ
    option="d"
    if(present(opt)) option=opt
    izero=0
    select case (option)
       case ("d")
          occ=0d0
          do jdg=1,nbasis
             if(abs(state(jdg)).lt.1d-15) cycle
             id = basis(jdg)
             !Extract bit information
             call magxtr(id,ksite,izero)
             do orb=1,6
                call cdag_i_c_j(ksite,ksitenew,INFO,orb,orb)
                if(INFO.eq.1) then
                   call magcon(ksitenew,ires,izero)
                   call binsrch(ires,idg,basis,nbasis)
                   call phase(ksite,i,j,phi)            	 
                  occ(orb)= occ(orb) + phi*conjg(state(idg))*state(jdg)
                endif
             enddo
          enddo
          do orb=1,6
             if(abs(imag(occ(orb))).gt.1d-6) then
                write(*,*) 'Complex occupancy is not allowed'
                stop
             endif
          enddo
          occupation=real(occ)
       case ("j")
          occ=0d0
          SO_evec=-HSO(1:6,1:6)
          call diagonalize(SO_evec,SO_eval)
          do jdg=1,nbasis
             if(abs(state(jdg)).lt.1d-15) cycle
             id = basis(jdg)
             !Extract bit information
             call magxtr(id,ksite,izero)
             do orb=1,6
                do i=1,6
                   do j=1,6
                      call cdag_i_c_j(ksite,ksitenew,INFO,i,j)
                      if(INFO.eq.1) then
                         call magcon(ksitenew,ires,izero)
                         call binsrch(ires,idg,basis,nbasis)
                         call phase(ksite,i,j,phi)            	 
                         occ(orb)= occ(orb) + phi*SO_evec(i,orb)*conjg(SO_evec(j,orb))*conjg(state(idg))*state(jdg)
                      endif
                   enddo
                enddo
             enddo
          enddo
          do orb=1,6
             if(abs(imag(occ(orb))).gt.1d-6) then
                write(*,*) 'Complex occupancy is not allowed'
                stop
             endif
          enddo
          occupation=real(occ)
    end select
  end subroutine single_particle_occupation



  real(kind=double) function charge_fluctuation(state,nbasis,basis)
    complex(kind=double),dimension(:),intent(in)::state
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    integer::i,j,jdg,idg,id,inew,ires,izero
    integer,dimension(nsite)::ksite
    real(kind=double)::avgN,avgN2,tmp
    avgN=0d0; avgN2=0d0
    do jdg=1,nbasis
       id = basis(jdg)
       !Extract bit information
       call magxtr(id,ksite,izero)
       tmp=sum(ksite(1:6))
       avgN=avgN+tmp*(abs(state(jdg))**2)
       avgN2=avgN2+(tmp**2)*(abs(state(jdg))**2)
    enddo
    charge_fluctuation=sqrt(avgN2-(avgN**2))
  end function charge_fluctuation

  subroutine jz_matrix(state,nbasis,basis,jzmat)
    implicit none
    complex(kind=double),dimension(:,:),intent(in)::state
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:,:),intent(out)::jzmat 
    if(j_init.ne.0) then
       call generate_j(nbasis,basis)
       j_init=0
    endif
    jzmat=matmul(transpose(conjg(state)),matmul(j_z,state))
  end subroutine jz_matrix

  real(kind=double) function j_moment(state,nbasis,basis)
    complex(kind=double),dimension(:,:),intent(in)::state
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:,:),allocatable::jmat
    complex(kind=double)::temp
    integer::i,njbasis
    if(j_init.ne.0) then
       call generate_j(nbasis,basis)
       j_init=0
    endif
    njbasis=size(state,2)
    allocate(jmat(njbasis,njbasis))
    jmat=matmul(transpose(conjg(state)),matmul(j2,state))
    do i=1,njbasis
       temp=temp+jmat(i,i)
    enddo
    if(imag(temp).gt.1d-7) then
       write(*,*) 'J matrix is not Hermitian'
       write(*,*) temp
       stop
    endif
    j_moment=real(temp)/njbasis
  end function j_moment

  real(kind=double) function s_moment(state,nbasis,basis)
    complex(kind=double),dimension(:),intent(in)::state
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:),allocatable::newstate
    complex(kind=double)::temp
    if(j_init.ne.0) then
       call generate_j(nbasis,basis)
       j_init=0
    endif
    if(allocated(newstate)) deallocate(newstate)
    allocate(newstate(nbasis))
    newstate=matmul(s2,state)
    temp=dot_product(state,newstate)
    s_moment=real(temp)
  end function s_moment

  real(kind=double) function l_moment(state,nbasis,basis)
    complex(kind=double),dimension(:),intent(in)::state
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:),allocatable::newstate
    complex(kind=double)::temp
    if(j_init.ne.0) then
       call generate_j(nbasis,basis)
       j_init=0
    endif
    if(allocated(newstate)) deallocate(newstate)
    allocate(newstate(nbasis))
    newstate=matmul(l2,state)
    temp=dot_product(state,newstate)
    l_moment=real(temp)
  end function l_moment


  real(kind=double) function gyro_ratio(states,nbasis,basis)
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:,:),intent(in)::states
    complex(kind=double),dimension(:,:),allocatable::tmp_mat,smat,lmat,jmat
    complex(kind=double)::g_tmp,moment,jval
    integer::i,j,nstate
    if(j_init.ne.0) then
       call generate_j(nbasis,basis)
       j_init=0
    endif
    nstate=size(states,2)
    allocate(smat(nstate,nstate),lmat(nstate,nstate),jmat(nstate,nstate),tmp_mat(nbasis,nstate))
    tmp_mat=matmul(j_z,states)
    jmat=matmul(conjg(transpose(states)),tmp_mat)
    tmp_mat=matmul(s_z,states)
    smat=matmul(conjg(transpose(states)),tmp_mat)
    tmp_mat=matmul(l_z,states)
    lmat=matmul(conjg(transpose(states)),tmp_mat)
    gyro_ratio=0d0
    !do i=1,nstate
    !   write(*,'(20(F6.3,2x))') (2.d0*smat(i,:)-lmat(i,:))
    !enddo
    !do i=1,nstate
    !   write(*,'(20(F6.3,2x))') smat(i,:)
    !enddo
    do i=1,nstate
       do j=1,nstate
          if(abs(jmat(i,j)).gt.1d-3) then
             gyro_ratio=0d0
             g_tmp=smat(i,j)/jmat(i,j)
             gyro_ratio=gyro_ratio+2.d0*real(g_tmp)
             g_tmp=lmat(i,j)/jmat(i,j)
             gyro_ratio=gyro_ratio-real(g_tmp)
             exit
          endif
       enddo
    enddo
  end function  gyro_ratio


  function siteJ_moment(state,nbasis,basis)
    implicit none
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:),intent(in)::state
    real(kind=double),dimension(3)::siteJ_moment
    complex(kind=double)::tmp
    if(j_init.ne.0) then
       call generate_j(nbasis,basis)
       j_init=0
    endif
    siteJ_moment=0d0
    tmp=dot_product(state,matmul(j_x,state))
    if(abs(imag(tmp)).gt.1d-5) write(*,*) 'Moment not real'
    siteJ_moment(1)=real(tmp)
    tmp=dot_product(state,matmul(j_y,state))
    if(abs(imag(tmp)).gt.1d-5) write(*,*) 'Moment not real'
    siteJ_moment(2)=real(tmp)
    tmp=dot_product(state,matmul(j_z,state))
    if(abs(imag(tmp)).gt.1d-5) write(*,*) 'Moment not real'
    siteJ_moment(3)=real(tmp)
  end function siteJ_moment


  subroutine generate_j(nbasis,basis)
    integer,intent(in)::nbasis
    integer,dimension(:),intent(in)::basis
    integer::i,j,k,l,jdg,idg,id,inew,inew2,ires,izero
    integer,dimension(nsite)::ksite,ksitenew,ksitenew2
    real(kind=double)::phi,phitot
    !Define single particle j operators
    jx=0d0; jy=0d0; jz=0d0
    jx(1,2)=0.5d0;       jx(2,1)=0.5d0;       jx(3,4)=0.5d0;       jx(3,5)=(0d0,-1.d0)
    jx(4,3)=0.5d0;       jx(4,6)=(0d0,-1.d0);  jx(5,3)=(0d0,1.d0); jx(5,6)=0.5d0
    jx(6,4)=(0d0,1.d0); jx(6,5)=0.5d0
    jz(1,1)=0.5d0;       jz(1,3)=(0d0,-1.d0);  jz(2,2)=-0.5d0;      jz(2,4)=(0d0,-1.d0)
    jz(3,1)=(0d0,1.d0); jz(3,3)=0.5d0;       jz(4,2)=(0d0,1.d0); jz(4,4)=-0.5d0
    jz(5,5)=0.5d0;       jz(6,6)=-0.5d0
    !Uncomment next two lines for j moments on site 2
    !jx(7:12,7:12)=jx(1:6,1:6); jx(1:6,1:6)=0d0
    !jz(7:12,7:12)=jz(1:6,1:6); jz(1:6,1:6)=0d0
    jy=(0d0,-1.d0)*(matmul(jz,jx)-matmul(jx,jz))
    !Use this for total J of the two site problem
    jxx=0d0; jyy=0d0; jzz=0d0
    jxx(1:6,1:6)=jx(1:6,1:6); jxx(7:12,7:12)=jx(1:6,1:6)
    jyy(1:6,1:6)=jy(1:6,1:6); jyy(7:12,7:12)=jy(1:6,1:6)
    jzz(1:6,1:6)=jz(1:6,1:6); jzz(7:12,7:12)=jz(1:6,1:6)
    !Define single particle s operators
    sx=0d0; sy=0d0; sz=0d0
    sx(1,2)=0.5d0;   sx(2,1)=0.5d0;   sx(3,4)=0.5d0;   sx(4,3)=0.5d0
    sx(5,6)=0.5d0;   sx(6,5)=0.5d0
    sz(1,1)=0.5d0;   sz(2,2)=-0.5d0;  sz(3,3)=0.5d0;   sz(4,4)=-0.5d0
    sz(5,5)=0.5d0;   sz(6,6)=-0.5d0
    sy=(0d0,-1.d0)*(matmul(sz,sx)-matmul(sx,sz))
    !Use this for total S of the two site problem
    sxx=0d0; syy=0d0; szz=0d0
    sxx(1:6,1:6)=sx(1:6,1:6);  sxx(7:12,7:12)=sx(1:6,1:6)
    syy(1:6,1:6)=sy(1:6,1:6);  syy(7:12,7:12)=sy(1:6,1:6)
    szz(1:6,1:6)=sz(1:6,1:6);  szz(7:12,7:12)=sz(1:6,1:6)
    !Define single particle l operators
    lx=0d0; ly=0d0; lz=0d0
    lx(3,5)=(0d0,-1.d0); lx(4,6)=(0d0,-1.d0); lx(5,3)=(0d0,1.d0);  lx(6,4)=(0d0,1.d0)
    ly(1,5)=(0d0,1.d0);  ly(2,6)=(0d0,1.d0);  ly(5,1)=(0d0,-1.d0); ly(6,2)=(0d0,-1.d0)
    lz=(0d0,-1.d0)*(matmul(lx,ly)-matmul(ly,lx))
    !Use this for total L of the two site problem
    lxx=0d0; lyy=0d0; lzz=0d0
    lxx(1:6,1:6)=lx(1:6,1:6);  lxx(7:12,7:12)=lx(1:6,1:6)
    lyy(1:6,1:6)=ly(1:6,1:6);  lyy(7:12,7:12)=ly(1:6,1:6)
    lzz(1:6,1:6)=lz(1:6,1:6);  lzz(7:12,7:12)=lz(1:6,1:6)  
    ! Generate J^2 operators
    allocate(j2(nbasis,nbasis),s2(nbasis,nbasis),l2(nbasis,nbasis))
    allocate(s_z(nbasis,nbasis),l_z(nbasis,nbasis),j_z(nbasis,nbasis),j_x(nbasis,nbasis),j_y(nbasis,nbasis))
    j2=0d0; s2=0d0; l2=0d0
    j_z=0d0; s_z=0d0; l_z=0d0
    j_x=0d0; j_y=0d0
    do jdg=1,nbasis
       id = basis(jdg)
       !Extract bit information
       call magxtr(id,ksite,izero)
       do i=1,12
          do j=1,12
             do k=1,12
                do l=1,12
                   phitot=1.d0
                   call cdag_i_c_j(ksite,ksitenew,inew,k,l) !cdag_up*c_dn
                   if(inew.eq.1) then
                      call phase(ksite,k,l,phi)
                      phitot=phitot*phi
                      call cdag_i_c_j(ksitenew,ksitenew2,inew2,i,j) !cdag_dn*c_up
                      if(inew2.eq.1) then
                         call magcon(ksitenew2,ires,izero)
                         call binsrch(ires,idg,basis,nbasis)
                         call phase(ksitenew,i,j,phi)
                         phitot=phitot*phi
                         !if(phitot.lt.0) write(*,*) i,j,k,l,ksite
                         !j2(idg,jdg)= j2(idg,jdg) + phitot*(jx(i,j)*jx(k,l) + jy(i,j)*jy(k,l) + jz(i,j)*jz(k,l))
                         !s2(idg,jdg)= s2(idg,jdg) + phitot*(sx(i,j)*sx(k,l) + sy(i,j)*sy(k,l) + sz(i,j)*sz(k,l))
                         !l2(idg,jdg)= l2(idg,jdg) + phitot*(lx(i,j)*lx(k,l) + ly(i,j)*ly(k,l) + lz(i,j)*lz(k,l))
                         j2(idg,jdg)= j2(idg,jdg) + phitot*(jxx(i,j)*jxx(k,l) + jyy(i,j)*jyy(k,l) + jzz(i,j)*jzz(k,l))
                         s2(idg,jdg)= s2(idg,jdg) + phitot*(sxx(i,j)*sxx(k,l) + syy(i,j)*syy(k,l) + szz(i,j)*szz(k,l))
                         l2(idg,jdg)= l2(idg,jdg) + phitot*(lxx(i,j)*lxx(k,l) + lyy(i,j)*lyy(k,l) + lzz(i,j)*lzz(k,l))
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       do i=1,12
          do j=1,12
             call cdag_i_c_j(ksite,ksitenew,inew,i,j) !cdag_up*c_dn
             if(inew.eq.1) then
                call phase(ksite,i,j,phi)
                call magcon(ksitenew,ires,izero)
                call binsrch(ires,idg,basis,nbasis)
                !l_z(idg,jdg)= l_z(idg,jdg) + phi*lz(i,j)
                !s_z(idg,jdg)= s_z(idg,jdg) + phi*sz(i,j)
                j_z(idg,jdg)= j_z(idg,jdg) + phi*jz(i,j)
                j_x(idg,jdg)= j_x(idg,jdg) + phi*jx(i,j)
                j_y(idg,jdg)= j_y(idg,jdg) + phi*jy(i,j)
                l_z(idg,jdg)= l_z(idg,jdg) + phi*lzz(i,j)
                s_z(idg,jdg)= s_z(idg,jdg) + phi*szz(i,j)
                !j_z(idg,jdg)= j_z(idg,jdg) + phi*jzz(i,j)                
             endif
          enddo
       enddo
    enddo
  end subroutine generate_j

  subroutine diagonalize_real(ham,eval)
    implicit none
    real(kind=double),dimension(:,:),intent(inout)::ham
    real(kind=double),dimension(:),intent(out)::eval
    integer::LWORK,ndim
    real(kind=double),dimension(:),allocatable::WORK
    integer::INFO
    ndim=size(ham,1)
    LWORK=3*ndim
    allocate(WORK(LWORK),stat=INFO)
    if(INFO.ne.0) then
       write(*,*) 'allocation of diagonalization arrays failed'
       stop
    endif
    call DSYEV('V', 'U', ndim, ham, ndim, eval, WORK, LWORK, INFO )
    if(INFO.ne.0) then
       write(*,*) 'diagonalization failed'
       stop
    endif
    deallocate(WORK)
    !write(*,*) 'Hamiltonian diagonalized'
  end subroutine diagonalize_real


  subroutine diagonalize_complex(ham,eval)
    implicit none
    complex(kind=double),dimension(:,:),intent(inout)::ham
    real(kind=double),dimension(:),intent(out)::eval
    integer::LWORK,ndim
    complex(kind=double),dimension(:),allocatable::WORK
    real(kind=double),dimension(:),allocatable::RWORK
    integer::INFO
    ndim=size(ham,1)
    LWORK=3*ndim
    allocate(WORK(LWORK),RWORK(3*ndim-2),stat=INFO)
    if(INFO.ne.0) then
       write(*,*) 'allocation of diagonalization arrays failed'
       stop
    endif
    call ZHEEV('V', 'U', ndim, ham, ndim, eval, WORK, LWORK, RWORK, INFO )
    if(INFO.ne.0) then
       write(*,*) 'diagonalization failed'
       stop
    endif
    deallocate(WORK,RWORK)
    !write(*,*) 'Hamiltonian diagonalized'
  end subroutine diagonalize_complex


  subroutine initialize(nbasis)
    !Initialize nmltplt,nbitsps,lbit and nto2,lmax and ibin
    implicit none
    integer,intent(in)::nbasis
    integer::i,iter,nbit
    nmltplt=2
    !nbitsps
    lmax=0
    do i=1,nsite
       nbit=0
       iter=1
5      iter=2*iter
       nbit=nbit+1
       if(iter.lt.(nmltplt(i))) goto 5
       nbitsps(i)=nbit
       lmax=lmax+nbitsps(i)
    enddo
    !lbit
    lbit(1)=lmax-nbitsps(1)
    do i=2,nsite
       lbit(i)=lbit(i-1)-nbitsps(i)
    enddo
    !nto2
    nto2(1)=1
    iter=1
6   iter=iter+1
    nto2(iter)=nto2(iter-1)*2
    if(nto2(iter).le.nbasis) goto 6
    ibin=iter
    !write(*,*) 'lmax,ibin',lmax,ibin
  end subroutine initialize

  subroutine magxtr(id,ksite,izero)
    implicit none
    integer,intent(in)::id,izero
    integer,dimension(:),intent(out)::ksite
   integer::j,k,ibeg
    do k=1,nsite
       ksite(k)=0
    enddo
    do j=1,nsite
       ibeg=lbit(j)
       call mvbits(id,ibeg,nbitsps(j),ksite(j),izero)
    enddo
  end subroutine magxtr

  subroutine magcon(ksitenew,ires,izero)
    implicit none
    integer,intent(in),dimension(:)::ksitenew
    integer,intent(in)::izero
    integer,intent(out)::ires
    integer::j,ibeg
    ires=0
    do j=1,nsite
       ibeg=lbit(j)
       call mvbits(ksitenew(j),izero,nbitsps(j),ires,lbit(j))
    enddo
  end subroutine magcon
      
  subroutine binsrch(id1,ifound,idet,islt)
    implicit none
    integer,intent(in)::id1,islt
    integer,intent(in),dimension(:)::idet
    integer,intent(out)::ifound
    integer::ibeg,itst
    ibeg=ibin-1
    itst=nto2(ibeg)
35  if(itst.gt.islt) go to 50
    if(id1.eq.idet(itst)) go to 25
    if(id1.gt.idet(itst)) go to 15
50  ibeg=ibeg-1
    if(ibeg.eq.0) go to 60
    itst=itst-nto2(ibeg)
    go to 35
15  ibeg=ibeg-1
    if(ibeg.eq.0) go to 60
    itst=itst+nto2(ibeg)
    go to 35
60  write(*,*)'error final state not found in 2nd diagram search'
    stop
25  ifound=itst
  end subroutine binsrch

  subroutine phase(ksite,i,j,phi)
    implicit none
    integer,intent(in)::i,j
    integer,intent(in),dimension(:)::ksite
    real(kind=double),intent(out)::phi
    integer::istart,iend,k
    phi=1.0d0
    if(i.eq.j) return
    if(i.gt.j) then
       istart=i
       iend=j+1	
    else
       istart=j
       iend=i+1
       phi=-1.0d0
    endif
    do k=istart,iend,-1
       if(ksite(k).eq.1) phi=phi*(-1.d0)
    enddo
  end subroutine phase

  subroutine cdag_i_c_j(ksite,ksitenew,inew,i,j)
    implicit none
    integer,intent(in),dimension(:)::ksite
    integer,intent(out),dimension(:)::ksitenew
    integer,intent(out)::inew
    integer,intent(in)::i,j
    inew=0
    ksitenew=ksite
    if(i.eq.j) then
       if(ksite(i).eq.1) inew=1
       return
    endif
    if((ksite(i).eq.1).or.(ksite(j).eq.0)) return
    inew=1
    ksitenew(i)=1
    ksitenew(j)=0
  end subroutine cdag_i_c_j

end module hamiltonian

 
