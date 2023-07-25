integer :: i, j, k
integer, allocatable :: grid2patch(:)
integer, allocatable :: grid2patch_start(:,:)
integer, allocatable :: grid2patch_count(:,:)

allocate(grid2patch(numpatch))
allocate(grid2patch_start(lon_points, lat_points))
allocate(grid2patch_count(lon_points, lat_points))

! calculate the number of patches of each gridpoint
grid2patch_count = 0
do k=1, numpatch
  i = ixy_patch(k)
  j = jxy_patch(k)
  grid2patch_count(i,j) = grid2patch_count(i,j) + 1
enddo

! calculate the start of patches of each gridpoint
k = 1
do j=1, lat_points
  do i=1, lon_points
    grid2patch_start(i,j) = k
    
    k = k + grid2patch_count(i,j)
  enddo
enddo

! put the patches of each gridpoint into
do k=1, numpatch
  i = ixy_patch(k)
  j = jxy_patch(k)

  grid2patch(grid2patch_start(i,j)) = k
  grid2patch_start(i,j) = grid2patch_start(i,j) + 1
enddo

! recover grid2patch_start
grid2patch_start = grid2patch_start - grid2patch_count

