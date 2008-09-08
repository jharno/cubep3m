! plot power spectra
! f95 driver.f90 plot_power.f90 `plplot-config --cflags --libs --with-f77`
subroutine plot_power(nc,nx,power,nm)
implicit none
integer, parameter :: PL_PARSE_FULL = 1
real*8, parameter :: PI = 3.1415926535897932384d0
integer, intent(in) :: nx,nc
real*8, dimension(nc,nx) :: power 
real*8 :: xr_min,xr_max,yr_min,yr_max     
character*80 :: fn
integer :: i,j,k
character*80 :: nm 

print *,minval(power(1,:))
print *,minval(power(2:nc,:))
print *,minval(power(2,:))
print *,'-----'
!convert power to log10
power=log10(power)

!Set up plot
call plsdev('png')
write(fn,'(a,a)') nm(1:len_trim(nm)-4),'.png'
call plsfnam(fn)
call plinit
call plfont(2)
call pladv(0)

! set the range of the plot from 0:1 in x and y in terms of total viewable area  
call plvpor(0.15d0, 0.85d0, 0.1d0, 0.9d0)

! set range of axis
xr_min=minval(floor(power(1,:)))
xr_max=maxval(ceiling(power(1,:)))
yr_min=minval(floor(power(2,:)))
yr_max=maxval(ceiling(power(2,:)))
print *,xr_min,xr_max,yr_min,yr_max
call plwind(xr_min, xr_max, yr_min, yr_max)
call plcol0(1)

!Try different axis and labelling styles.
call plbox('bclnstgh', 0.0d0, 0, 'bclnstvg', 0.0d0, 0)

!Plot spectra
do i=2,nc
  call plcol0(i+1)
  call plline(nx,power(1,:),power(i,:))
enddo

!Put labels on.
call plcol0(3)
call plmtex('b', 3.2d0, 0.5d0, 0.5d0, 'k (h/Mpc)')
call plmtex('t', 2.0d0, 0.5d0, 0.5d0,'Power spectrum')
call plcol0(2)
call plmtex('l', 5.0d0, 0.5d0, 0.5d0, 'Delta*Delta')
call plend

end subroutine plot_power
