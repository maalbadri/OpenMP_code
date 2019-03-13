program area_mandelbrot

  use omp_lib
  implicit none
!
  integer, parameter :: sp = kind(1.0)   
  integer, parameter :: dp = kind(1.0d0)
  integer :: i, j, iter, numoutside, id
  integer, parameter :: npoints = 2000, maxiter = 2000
  real (kind=dp) :: area, error, wtime
  complex (kind=dp) :: c , z
  integer :: threads, max_threads

max_threads = omp_get_max_threads ( )
!
! Calculate area of mandelbrot set
!
!    Outer loops runs over npoints, initialize z=c
!
!    Inner loop has the iteration z=z*z+c, and threshold test
!

write(*,*) ' __       __                            __            __  __                              __     ' 
write(*,*) '/  \     /  |                          /  |          /  |/  |                            /  |    '
write(*,*) '$$  \   /$$ |  ______   _______    ____$$ |  ______  $$ |$$ |____    ______    ______   _$$ |_   '
write(*,*) '$$$  \ /$$$ | /      \ /       \  /    $$ | /      \ $$ |$$      \  /      \  /      \ / $$   |  '
write(*,*) '$$$$  /$$$$ | $$$$$$  |$$$$$$$  |/$$$$$$$ |/$$$$$$  |$$ |$$$$$$$  |/$$$$$$  |/$$$$$$  |$$$$$$/   '
write(*,*) '$$ $$ $$/$$ | /    $$ |$$ |  $$ |$$ |  $$ |$$    $$ |$$ |$$ |  $$ |$$ |  $$/ $$ |  $$ |  $$ | __ '
write(*,*) '$$ |$$$/ $$ |/$$$$$$$ |$$ |  $$ |$$ \__$$ |$$$$$$$$/ $$ |$$ |__$$ |$$ |      $$ \__$$ |  $$ |/  |'
write(*,*) '$$ | $/  $$ |$$    $$ |$$ |  $$ |$$    $$ |$$       |$$ |$$    $$/ $$ |      $$    $$/   $$  $$/ '
write(*,*) '$$/      $$/  $$$$$$$/ $$/   $$/  $$$$$$$/  $$$$$$$/ $$/ $$$$$$$/  $$/        $$$$$$/     $$$$/  '

write(*,*) ' '
write(*,*) 'Please specify # threads:'
write(*,*) ' '
read(*,*) threads
if (threads > max_threads) then
  write(*,*) ' '
  write(*,*) 'Too many threads you silly goose > . > '
  write(*,*) ' '
  threads = max_threads
  write(*,*) 'Running with your maximum number of threads:', max_threads
end if

  numoutside = 0 

call omp_set_num_threads(threads)
  wtime = omp_get_wtime()
!$omp parallel do &
!$omp shared (threads) &
!$omp private (z,c,iter) & 
!$omp reduction(+:numoutside) default(none)

  do i = 0,npoints-1 
     do j= 0,npoints-1 
        c = cmplx(-2.0+(2.5*i)/npoints + 1.0d-07,(1.125*j)/npoints + 1.0d-07)
        z = c
        iter = 0 
        do while (iter < maxiter) 
           z = z*z + c 
           iter = iter + 1
           if (real(z)*real(z)+aimag(z)*aimag(z) > 4) then
              numoutside = numoutside + 1 
              exit
           endif
        end do 
     end do
  end do
!!$omp end do
!!$omp end parallel


!  Finish up by measuring the elapsed time.
!
  wtime = omp_get_wtime ( ) - wtime

  write ( *, * ) ' ' 
  write ( *, * ) '  Elapsed wall clock time = ', wtime
  write ( *, * ) ' ' 
!
!  Terminate.


!
! Output results
!
  area = 2.0*2.5*1.125 * real(npoints*npoints-numoutside)/real(npoints*npoints)
  error = area/real(npoints)
  print *, "Area of Mandelbrot set = ",area," +/- ",error
!
  stop
end program area_mandelbrot
