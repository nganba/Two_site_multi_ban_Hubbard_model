! two_site_RXS.f90
! 
! Generates the two site RXS scattering amplitude matrix at a given energy
!
! Author: Oinam Nganba Meetei
! Date:   10/11/13
!
!******************************************************************************************

module RXS
  use atomic_basis
  use atomic_hamiltonian
  use basis_states
  use hamiltonian
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer,parameter:: nsite=6
  integer:: lmax,GSibin,Exibin
  integer,save::init=-1
  integer,dimension(nsite)::lbit,nbitsps,GSnto2,Exnto2,nmltplt
  public:: RXS_scattering_amplitude

contains
  subroutine RXS_scattering_amplitude(U,JH,Lambda,GS_Nel,L2edge,L3edge)
    real(kind=double),intent(in):: U,JH,Lambda
    integer,intent(in)::GS_Nel
    complex(kind=double),dimension(3,3),intent(out)::L2edge,L3edge
    complex(kind=double),dimension(:,:),allocatable::GS_evec,Ex_evec,den_mat,Lmat,jstate,jzmat
    real(kind=double),dimension(:),allocatable::GS_eval,Ex_eval,jz
    integer,dimension(1000)::GS_basis,Ex_basis,GSatomic_basis
    complex(kind=double),dimension(6,6)::p_states,ir_mat,jr_mat,Hso
    real(kind=double),dimension(6)::soc_eval
    real(kind=double),dimension(3,6,6)::rvec
    integer,dimension(2*nsite)::ksitetot1,ksitetot2
    integer,dimension(nsite)::ksite1,ksite2
    integer::id,idg,jdg,jdg1,jdg2,ires,izero,inew,Ex_Nel,GSatomic_Nel,GS_nbasis,Ex_nbasis,GSatomic_nbasis
    integer::i,j,k,l,GSdeg,Exdeg,iGS,iEx,iPol,jPol,iP,nbit,ibeg
    complex(kind=double)::tmp1,tmp2
    real(kind=double)::phi
    
    !GSatomic_Nel=1; Ex_Nel=2

    !Generate states in GS manifold
    call generate_basis(GS_Nel,GS_nbasis,GS_basis)
    allocate(GS_evec(GS_nbasis,GS_nbasis),GS_eval(GS_nbasis))
    call generate_hamiltonian(U,JH,Lambda,GS_Nel,GS_nbasis,GS_basis,GS_evec)
    call diagonalize(GS_evec,GS_eval)
    GSdeg=0
    do i=1,GS_nbasis
       if(abs(GS_eval(i)-GS_eval(1)).lt.1d-5) then
          GSdeg=GSdeg+1
       else
          exit
       endif
    enddo
    if(GSdeg.gt.1) then
       allocate(jstate(GS_nbasis,GSdeg),jzmat(GSdeg,GSdeg),jz(GSdeg))
       jstate=GS_evec(:,1:GSdeg)
       call jz_matrix(jstate,GS_nbasis,GS_basis,jzmat)
       call diagonalize(jzmat,jz)
       do i=1,GSdeg
          GS_evec(:,i)=0d0
          do k=1,GSdeg
             GS_evec(:,i)=GS_evec(:,i)+jzmat(k,i)*jstate(:,k)
          enddo
       enddo
       deallocate(jstate,jzmat,jz)
    endif

    !Generate excited state manifold
    Ex_Nel=(GS_Nel/2)+1
    write(*,*) 'Number of el in excited state', Ex_Nel
    call generate_atomic_basis(Ex_Nel,Ex_nbasis,Ex_basis)
    allocate(Ex_evec(Ex_nbasis,Ex_nbasis),Ex_eval(Ex_nbasis))
    call generate_atomic_H(U,JH,Lambda,Ex_Nel,Ex_nbasis,Ex_basis,Ex_evec)
    call diagonalize(Ex_evec,Ex_eval)
    Exdeg=0
    do i=1,Ex_nbasis
       if(abs(Ex_eval(i)-Ex_eval(1)).lt.1d-5) then
          Exdeg=Exdeg+1
       else
          exit
       endif
    enddo

    !Generate GS atomic basis states
    GSatomic_Nel=GS_Nel/2
    call generate_atomic_basis(GSatomic_Nel,GSatomic_nbasis,GSatomic_basis)

    !Initialize
    call initialize(GSatomic_nbasis,Ex_nbasis)


    !Calculate reduced density matrix for site 1
    izero=0
    allocate(den_mat(GSatomic_nbasis,GSatomic_nbasis))
    den_mat=0d0
    iGS=GSdeg
    write(*,*) 'GSdeg',GSdeg
    nbit=1
    !do iGS=1,GSdeg
    do i=1,GS_nbasis
       id=GS_basis(i)
       ksitetot1=0
       ibeg=12
       do k=1,12
          ibeg=ibeg-1
          call mvbits(id,ibeg,nbit,ksitetot1(k),izero)
       enddo
       !write(*,'(12(I3,1x))') ksitetot1
       if(sum(ksitetot1(1:6)).ne.GSatomic_Nel) cycle
       inner: do j=1,GS_nbasis
          id=GS_basis(j)
          ksitetot2=0
          ibeg=12
          do k=1,12
             ibeg=ibeg-1
             call mvbits(id,ibeg,nbit,ksitetot2(k),izero)
          enddo
          if(sum(ksitetot2(1:6)).ne.GSatomic_Nel) cycle inner
          do k=7,12
             if(ksitetot1(k).ne.ksitetot2(k)) cycle inner
          enddo
          !write(*,'(12(I3,2x))') i,ksitetot1
          !write(*,'(12(I3,2x))') j,ksitetot2
          !write(*,*)
          ksite1=ksitetot1(1:6)
          ksite2=ksitetot2(1:6)
          call magcon(ksite1,ires,izero)
          call GS_binsrch(ires,idg,GSatomic_basis,GSatomic_nbasis) 
          call magcon(ksite2,ires,izero)
          call GS_binsrch(ires,jdg,GSatomic_basis,GSatomic_nbasis)
          den_mat(idg,jdg)=den_mat(idg,jdg) + GS_evec(i,iGS)*conjg(GS_evec(j,iGS))
       enddo inner
    enddo
    !enddo
    !den_mat=den_mat/GSdeg
    do i=1,GSatomic_nbasis
       do j=i,GSatomic_nbasis
          if(abs(den_mat(i,j)-conjg(den_mat(j,i))).gt.1d-5) then
             write(*,*) 'Not hermitian density matrix'
             stop
          endif
       enddo
    enddo
    !write(*,*) 'den_mat'
    !do i=1,GSatomic_nbasis
    !   write(*,'(20(F5.2,1x,F5.2,2x))') den_mat(i,1:10)
    !enddo
    !tmp1=0d0
    !do i=1,GSatomic_nbasis
    !   tmp1=tmp1+den_mat(i,i)
    !enddo
    !write(*,*) 'Trace den_mat',tmp1
    
    !write(*,*) 'Generated reduced density matrix'

    !Calculate dipole excitation matrix elements
    !Initialize matrix elements of dipole operator
    ! Note: Basis used is in the order {dyz_up, dyz_dn, dzx_up, dzx_dn, dxy_up, dxy_dn}
    !       p-orbitals are ordered as {x_up, x_dn, y_up, y_dn, z_up, z_dn}
    rvec=0d0
    !Polarization along x
    rvec(1,3,5)=1.d0; rvec(1,4,6)=1.d0; rvec(1,5,3)=1.d0; rvec(1,6,4)=1.d0
    !Polarization along y
    rvec(2,1,5)=1.d0; rvec(2,2,6)=1.d0; rvec(2,5,1)=1.d0; rvec(2,6,2)=1.d0
    !Polarization along z
    rvec(3,1,3)=1.d0; rvec(3,2,4)=1.d0; rvec(3,3,1)=1.d0; rvec(3,4,2)=1.d0

    !Generate p-states with a hole in either P1/2 and P3/2 states
    Hso=0d0; p_states=0d0
    Hso(1,3)=-(0.d0,1.d0);  Hso(1,6)=1.d0
    Hso(2,4)=(0.d0,1.d0);   Hso(2,5)=-1.d0
    Hso(3,1)=(0.d0,1.d0);   Hso(3,6)=-(0.d0,1.d0)
    Hso(4,2)=-(0.d0,1.d0);  Hso(4,5)=-(0.d0,1.d0)
    Hso(5,2)=-1.d0;         Hso(5,4)=(0.d0,1.d0)
    Hso(6,1)=1.d0;          Hso(6,3)=(0.d0,1.d0)
    call diagonalize(Hso,soc_eval)
    p_states=0d0
    do i=1,6
       do j=1,6
          p_states(j,i)=Hso(j,i)
       enddo
    enddo

    !Calculate RXS scattering amplitude matrix
    L2edge=0d0; L3edge=0d0
    allocate(Lmat(GSatomic_nbasis,GSatomic_nbasis))
    !Loop over excited states
    do iEx=3,6
       !L2 edge : sum over p-states
       do iP=1,2
          !Loop over polarizations
          do iPol=1,3
             do jPol=1,3
                ir_mat=rvec(iPol,:,:)
                jr_mat=rvec(jPol,:,:)
                Lmat=0d0
                do jdg1=1,GSatomic_nbasis
                   id=GSatomic_basis(jdg1)
                   call magxtr(id,ksite1,izero)
                   tmp1=0d0
                   do i=1,6
                      call cdag_i(ksite1,ksite2,inew,i)
                      if(inew.eq.1) then
                         call magcon(ksite2,ires,izero)
                         call Ex_binsrch(ires,idg,Ex_basis,Ex_nbasis)
                         call phase(ksite1,i,phi)  
                         do j=1,6
                            tmp1=tmp1+conjg(Ex_evec(idg,iEx))*p_states(j,iP)*phi*ir_mat(i,j)
                         enddo
                      endif
                   enddo
                   do jdg2=1,GSatomic_nbasis
                      id=GSatomic_basis(jdg2)
                      call magxtr(id,ksite1,izero)
                      tmp2=0d0
                      do i=1,6
                         call cdag_i(ksite1,ksite2,inew,i)
                         if(inew.eq.1) then
                            call magcon(ksite2,ires,izero)
                            call Ex_binsrch(ires,idg,Ex_basis,Ex_nbasis)
                            call phase(ksite1,i,phi)  
                            do j=1,6
                               tmp2=tmp2+conjg(Ex_evec(idg,iEx))*p_states(j,iP)*phi*jr_mat(i,j)
                            enddo
                         endif
                      enddo
                      Lmat(jdg1,jdg2)=conjg(tmp1)*tmp2
                   enddo
                enddo
                Lmat=matmul(den_mat,Lmat)
                do i=1,GSatomic_nbasis
                   L2edge(iPol,jPol)=L2edge(iPol,jPol)+Lmat(i,i)
                enddo
             enddo
          enddo
       enddo
       !L3 edge : sum over p-states
       do iP=3,6
          !Loop over polarizations
          do iPol=1,3
             do jPol=1,3
                ir_mat=rvec(iPol,:,:)
                jr_mat=rvec(jPol,:,:)
                Lmat=0d0
                do jdg1=1,GSatomic_nbasis
                   id=GSatomic_basis(jdg1)
                   call magxtr(id,ksite1,izero)
                   tmp1=0d0
                   do i=1,6
                      call cdag_i(ksite1,ksite2,inew,i)
                      if(inew.eq.1) then
                         call magcon(ksite2,ires,izero)
                         call Ex_binsrch(ires,idg,Ex_basis,Ex_nbasis)
                         call phase(ksite1,i,phi)  
                         do j=1,6
                            tmp1=tmp1+conjg(Ex_evec(idg,iEx))*p_states(j,iP)*phi*ir_mat(i,j)
                         enddo
                      endif
                   enddo
                   do jdg2=1,GSatomic_nbasis
                      id=GSatomic_basis(jdg2)
                      call magxtr(id,ksite1,izero)
                      tmp2=0d0
                      do i=1,6
                         call cdag_i(ksite1,ksite2,inew,i)
                         if(inew.eq.1) then
                            call magcon(ksite2,ires,izero)
                            call Ex_binsrch(ires,idg,Ex_basis,Ex_nbasis)
                            call phase(ksite1,i,phi)  
                            do j=1,6
                               tmp2=tmp2+conjg(Ex_evec(idg,iEx))*p_states(j,iP)*phi*jr_mat(i,j)
                            enddo
                         endif
                      enddo
                      Lmat(jdg1,jdg2)=conjg(tmp1)*tmp2
                   enddo
                enddo
                Lmat=matmul(den_mat,Lmat)
                do i=1,GSatomic_nbasis
                   L3edge(iPol,jPol)=L3edge(iPol,jPol)+Lmat(i,i)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine RXS_scattering_amplitude



  subroutine initialize(GSnbasis,Exnbasis)
    !Initialize nmltplt,nbitsps,lbit and nto2,lmax and ibin
    implicit none
    integer,intent(in)::GSnbasis,Exnbasis
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
    !nto2 GS
    GSnto2(1)=1
    iter=1
6   iter=iter+1
    GSnto2(iter)=GSnto2(iter-1)*2
    if(GSnto2(iter).le.GSnbasis) goto 6
    GSibin=iter
    !nto2 Ex
    Exnto2(1)=1
    iter=1
7   iter=iter+1
    Exnto2(iter)=Exnto2(iter-1)*2
    if(Exnto2(iter).le.Exnbasis) goto 7
    Exibin=iter
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
      
  subroutine GS_binsrch(id1,ifound,idet,islt)
    implicit none
    integer,intent(in)::id1,islt
    integer,intent(in),dimension(:)::idet
    integer,intent(out)::ifound
    integer::ibeg,itst
    ibeg=GSibin-1
    itst=GSnto2(ibeg)
35  if(itst.gt.islt) go to 50
    if(id1.eq.idet(itst)) go to 25
    if(id1.gt.idet(itst)) go to 15
50  ibeg=ibeg-1
    if(ibeg.eq.0) go to 60
    itst=itst-GSnto2(ibeg)
    go to 35
15  ibeg=ibeg-1
    if(ibeg.eq.0) go to 60
    itst=itst+GSnto2(ibeg)
    go to 35
60  write(*,*)'error final state not found in 2nd diagram search'
    stop
25  ifound=itst
  end subroutine GS_binsrch

  subroutine Ex_binsrch(id1,ifound,idet,islt)
    implicit none
    integer,intent(in)::id1,islt
    integer,intent(in),dimension(:)::idet
    integer,intent(out)::ifound
    integer::ibeg,itst
    ibeg=Exibin-1
    itst=Exnto2(ibeg)
35  if(itst.gt.islt) go to 50
    if(id1.eq.idet(itst)) go to 25
    if(id1.gt.idet(itst)) go to 15
50  ibeg=ibeg-1
    if(ibeg.eq.0) go to 60
    itst=itst-Exnto2(ibeg)
    go to 35
15  ibeg=ibeg-1
    if(ibeg.eq.0) go to 60
    itst=itst+Exnto2(ibeg)
    go to 35
60  write(*,*)'error final state not found in 2nd diagram search'
    stop
25  ifound=itst
  end subroutine Ex_binsrch
  

  subroutine phase(ksite,i,phi)
    implicit none
    integer,intent(in)::i
    integer,intent(in),dimension(:)::ksite
    real(kind=double),intent(out)::phi
    integer::istart,iend,k
    phi=1.d0
    if(i.eq.nsite) return
    k=sum(ksite(i+1:nsite))
    phi=(-1.d0)**k
  end subroutine phase

  subroutine cdag_i(ksite,ksitenew,inew,i)
    implicit none
    integer,intent(in),dimension(:)::ksite
    integer,intent(out),dimension(:)::ksitenew
    integer,intent(out)::inew
    integer,intent(in)::i
    inew=0
    ksitenew=ksite
    if((ksite(i).eq.1)) return
    inew=1
    ksitenew(i)=1
  end subroutine cdag_i

  subroutine c_i(ksite,ksitenew,inew,i)
    implicit none
    integer,intent(in),dimension(:)::ksite
    integer,intent(out),dimension(:)::ksitenew
    integer,intent(out)::inew
    integer,intent(in)::i
    inew=0
    ksitenew=ksite
    if((ksite(i).eq.0)) return
    inew=1
    ksitenew(i)=0
  end subroutine c_i 

end module RXS

 
