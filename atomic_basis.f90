! two_site_basis.f90      
!
! This program generates the two site basis states of a 3 band fermionic Huubard model for a given number 
! occupation number. 
! It is based on basis_ferion.f
!
! Author:    	Oinam Nganba Meetei
! Date: 	08/02/13 
!***************************************************************************************

module atomic_basis
  implicit none
  private
  integer,parameter::single = selected_real_kind(p=6,r=37)
  integer,parameter::double = selected_real_kind(p=13,r=200)
  integer,parameter:: maxel=1000, kci=100, mtmp=100, nsite=6
  integer,dimension(nsite)::nmltplt,nbitsps,ksite,lbit
  integer,dimension(maxel)::idgseq

  public:: generate_atomic_basis

contains

  subroutine generate_atomic_basis(nel,nbasis,basis)
    integer,intent(in)::nel
    integer,intent(out)::nbasis
    integer,dimension(:),intent(out)::basis
    integer:: lmax,nbit,iter,izero,imin,imax,id,idiag,k,ntot,j,ibeg,i
    ! nel               number of electrons
    ! nmltplt(i)	number of possible states on i th site
    ! nbitsps(i)	number of bits reserved for i th site
    ! idgseq(i) 	number representation of i th basis state
    ! ksite(i) 	        state of i th site in a given basis state
    ! lbit(i)	        bit position from which information about i th site starts
    ! nsite             number of orbitals (count up and down as separate)
    nmltplt=2  ! 2 options per orbital: empty or occupied
    lmax=0 ! lmax = total number of bits for all sites
    do i=1,nsite
       nbit=0
       iter=1
10     iter=2*iter
       nbit=nbit+1
       if(iter.lt.(nmltplt(i))) goto 10
       nbitsps(i)=nbit
       lmax=lmax+nbitsps(i)
    enddo    
    lbit(1)=lmax-nbitsps(1)
    do i=2,nsite
       lbit(i)=lbit(i-1)-nbitsps(i)
    enddo
    izero=0
    imin=0
    imax=(2**lmax) -1 !imax=largest integer to be sampled
    !write(*,*)'imin,imax',imin,imax
    idiag=0
    out: do id=imin,imax
       !if(id/10000000*10000000.eq.id) write(*,*)'i,idiag',id,idiag
       ksite=0
       ntot=0
       do j=1,nsite
          ibeg=lbit(j)
          call mvbits(id,lbit(j),nbitsps(j),ksite(j),izero)
          ! mvbits breaks id into bits and extract bit information starting from bit position
          ! lbit(j) to lbit(j)+nbitsps(j) and puts them to ksite(j)	  
          if(ksite(j).gt.(nmltplt(j)-1)) cycle out
	  if(ksite(j).eq.1) ntot=ntot+1
       enddo
       if(ntot.ne.nel) cycle out
       !write(*,*) id
       !write(*,'(20(I3,2x))') ksite
       idiag=idiag+1
       idgseq(idiag)=id
    enddo out
    basis=idgseq
    nbasis=idiag
  end subroutine generate_atomic_basis

end module atomic_basis


      

      
