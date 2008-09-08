implicit none

integer, parameter :: max_halos=128**3
real(4), dimension(6,max_halos) :: halo_list
integer(4), dimension(max_halos) :: isorthalo
real(4), dimension(max_halos) :: mass
integer::k,nh_sum
common / rarr / halo_list,mass
common / iarr / isorthalo

call random_seed
nh_sum=max_halos !400000
call random_number(halo_list(:,:nh_sum))
mass(:nh_sum)=halo_list(4,:nh_sum)
isorthalo(:nh_sum)=(/ (k,k=1,nh_sum) /)
call indexedsort(nh_sum,mass,isorthalo)
halo_list(:,:nh_sum)=halo_list(:,isorthalo(:nh_sum))
do k=nh_sum,nh_sum-10,-1
  print *,halo_list(:,k)
enddo
print *,halo_list(:,1)
end    
