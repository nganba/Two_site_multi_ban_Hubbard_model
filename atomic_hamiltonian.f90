! atomic_hamiltonian.f90
! 
! Generates the atomic Hamiltonian including the effects of U, J_H and spin orbit 
! in a 3 band model. It is based on hmatgen_fermion.f90
!
! It also has subroutines to diagonalize the Hamiltonian (copied from hamiltonian.f90)
!
! Author: Oinam Nganba Meetei
! Date:   08/02/13
!
!******************************************************************************************

module atomic_hamiltonian
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer,parameter:: nsite=6
  integer:: lmax,ibin
  integer,save::j_init=-1
  integer,dimension(nsite)::lbit,nbitsps,nto2,nmltplt
  complex(kind=double),dimension(6,6),save::Hso
  public:: generate_atomic_H

contains
  subroutine generate_atomic_H(U,JH,lambda,Nel,nbasis,basis,hmat)
    real(kind=double),intent(in):: U,JH,lambda
    integer,intent(in)::Nel,nbasis
    integer,dimension(:),intent(in)::basis
    complex(kind=double),dimension(:,:),intent(out)::hmat
    integer,dimension(nsite)::ksite,ksitenew,ksitenew2
    integer::id,jdg,i,j,k,idg,ires,izero,nbnd,nbonds,inew,inew2
    real(kind=double)::phi,phitot
    call initialize(nbasis)
    if(size(hmat,1).ne.nbasis) then
       write(*,*) 'Hamiltonian array needs to have the same size a basis states'
       stop
    endif
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
       !Tetragonal distortion
       !hmat(jdg,jdg)=hmat(jdg,jdg) + 0.5d0*(ksite(5)+ksite(6))
       !Interaction
       !Intra-orbital
       do i=1,6,2
          j=i+1
          hmat(jdg,jdg)= hmat(jdg,jdg) + U*ksite(i)*ksite(j)
       enddo
       !Inter-orbital (same spin)
       do i=1,4
          do j=i+2,6,2
             hmat(jdg,jdg)= hmat(jdg,jdg) + (U-3.d0*JH)*ksite(i)*ksite(j)
          enddo
       enddo
       !Inter-orbital (opposite spin)
       do i=1,4
          k=i+3
          if(mod(i,2).eq.0) k=i+1
          do j=k,6,2
             hmat(jdg,jdg)= hmat(jdg,jdg) + (U-2.d0*JH)*ksite(i)*ksite(j)
          enddo
       enddo
       !Spin-flip term
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
       !Pair Hopping term
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
    enddo
    !do i=1,nbasis
    !   write(*,'(10(F6.3,1x,F6.3,3x))') hmat(i,:)
    !enddo
    !write(*,*)
  end subroutine generate_atomic_H

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
    phi=1.d0
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
    ksitenew(j)=0
    ksitenew(i)=1
  end subroutine cdag_i_c_j

end module atomic_hamiltonian



 
    
    
 

    

  
	
 
