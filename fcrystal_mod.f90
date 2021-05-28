subroutine sp_seeds(n_new_seeds, new_centers,  new_t0s, &
			        G, r_min, n_seeds, centers, t0s)
  ! This subroutine checks if seeds are born, or killed by other
  ! previously born particles.
  !
  ! n_new_seeds
  ! new_centers - array of dimension (n_new_seeds, 3)
  ! new_t0s - birth time, must be ordered! Array of dim (n_new_seeds)
  ! G - growth speed. Scalar.
  ! r_min - distance from the edge of the spherulite, where already no new 
  !         spherulites are born
  !         This parameter helps to decrease the number of very small spherulites.
  !
  ! n_seeds
  ! centers (n_new_seeds + n_old_seeds, 3)
  ! t0s (n_new_seeds + n_old_seeds)
  !
  ! returns positions and birthtime of spontaneously born spheroulites
  !
  
  integer, intent(in) :: n_new_seeds
  real, dimension(n_new_seeds), intent(in) :: new_t0s
  real, dimension(n_new_seeds,3), intent(in) :: new_centers
  
  real, intent(in) :: G
  real, intent(in) :: r_min
  
  integer, intent(out) :: n_seeds
  real, dimension(n_new_seeds), intent(out) :: t0s
  real, dimension(n_new_seeds,3), intent(out) :: centers

  integer :: it
  logical :: dead
  real, dimension(n_new_seeds) :: dd
  
  n_seeds = 1
  centers(1,:) = new_centers(1,:)
  t0s(1) = new_t0s(1)

  do it = 2, n_new_seeds
    dead = .false.
    call arr_distance(new_centers(it,:), centers, n_seeds, dd(1:n_seeds))
    if (any(dd(1:n_seeds) < r_min + G*(new_t0s(it)-t0s))) then
      dead = .true.
    end if

    
     if (.not. dead) then
        n_seeds = n_seeds + 1
        centers(n_seeds,:) = new_centers(it,:)
        t0s(n_seeds) = new_t0s(it)
     end if
  end do
end subroutine sp_seeds

subroutine sp_seeds_aniso(n_new_seeds, new_centers,  new_t0s, &
			        T_0, Cr, params, r_min, n_seeds, centers, t0s)
  ! same as sp_seeds but for the anisotherm case
  ! n_new_seeds
  ! new_centers - array of dimension (n_new_seeds, 3)
  ! new_t0s - birth time, must be ordered! Array of dim (n_new_seeds)
  ! T_0 - starting temperature of crystallization [K]
  ! Cr: cooling rate [K/s]
  ! params: fitted parameters of the Hoffmann-Lauritzen equation
  ! r_min - distance from the edge of the spherulite, where already no new 
  !         spherulites are born
  !         This parameter helps to decrease the number of very small spherulites.
  !
  ! n_seeds
  ! centers (n_new_seeds + n_old_seeds, 3)
  ! t0s (n_new_seeds + n_old_seeds)
  !
  ! returns positions and birthtime of spontaneously born spheroulites
  
  integer, intent(in) :: n_new_seeds
  real, dimension(n_new_seeds), intent(in) :: new_t0s
  real, dimension(n_new_seeds,3), intent(in) :: new_centers
  
  real, intent(in) :: T_0, Cr
  real, intent(in) :: r_min
  
  real, dimension(2), intent(in) :: params
  
  integer, intent(out) :: n_seeds
  real, dimension(n_new_seeds), intent(out) :: t0s
  real, dimension(n_new_seeds,3), intent(out) :: centers

  integer :: it
  logical :: dead
  real, dimension(n_new_seeds) :: dd
  real :: temp_dd
  
  n_seeds = 1
  centers(1,:) = new_centers(1,:)
  t0s(1) = new_t0s(1)

  do it = 2, n_new_seeds
    dead = .false.
    call arr_distance(new_centers(it,:), centers, n_seeds, dd(1:n_seeds))
    call dist_calc(new_t0s(it)-t0s(it), t0s(it), T_0, Cr, params, temp_dd) !Taylor-sor
    if (any(dd(1:n_seeds) < r_min + temp_dd)) then
      dead = .true.
    end if


    
     if (.not. dead) then
        n_seeds = n_seeds + 1
        centers(n_seeds,:) = new_centers(it,:)
        t0s(n_seeds) = new_t0s(it)
     end if
  end do
end subroutine sp_seeds_aniso


subroutine sr_cdf(n_seeds, centers, t0s, LL, G, n_rnd_pos, t_conv, vols)
  ! calculates the cdf of a single realization, given the positions and 
  ! starting time of the seeds. 
  ! The subroutine assumes that all the seeds are born. Thus please run 
  ! sp_seeds subroutine before this one, if calculating spontaneous crystal growth.
  !
  ! n_seeds: number of spherulite seeds to grow
  ! centers: coordinates of the centers of spherulites (n_seeds,3)
  ! t0s: birthtime of the spherulites (n_seeds)
  ! LL: length of computational cube
  ! G: growing speed, scalar
  ! n_rnd_pos: number of random positions to check for crystallization time
  ! t_conv: to be returned, time of conversion at the positions (n_rnd_pos)
  ! vols: approximate volumes of the spherulites
  
  integer, intent(in) :: n_seeds, n_rnd_pos
  real, dimension(n_seeds, 3), intent(in) :: centers
  real, dimension(n_seeds), intent(in) :: t0s
  real, intent(in) :: LL, G
  real, dimension(n_rnd_pos), intent(out) :: t_conv
  real, dimension(n_seeds), intent(out) :: vols
  
  real, dimension(n_rnd_pos, 3) :: test_points
  real, dimension(n_seeds) :: t_conv_arr
  
  integer :: i, mloc

  real, dimension(n_seeds) :: dd_arr

  call init_random_seed()
  call random_number(test_points)
  test_points = test_points*LL
  vols = 0.0
  
  do i = 1, n_rnd_pos
    call arr_distance(test_points(i,:), centers, n_seeds, dd_arr)
    t_conv_arr = dd_arr/G+t0s
    mloc = minloc(t_conv_arr,1)
    t_conv(i) = t_conv_arr(mloc)
    vols(mloc) = vols(mloc)+1.0


 end do

 vols = vols/n_rnd_pos*LL**3
end subroutine sr_cdf


subroutine sr_cdf2(n_seeds, centers, t0s, LL, G, n_rnd_pos, t_conv, vols)
  ! calculates the cdf of a single realization, given the positions and 
  ! starting time of the seeds. 
  ! The subroutine assumes that all the seeds are born. Thus please run 
  ! sp_seeds subroutine before this one, if calculating spontaneous crystal growth.
  !
  ! n_seeds: number of spherulite seeds to grow
  ! centers: coordinates of the centers of spherulites (n_seeds,3)
  ! t0s: birthtime of the spherulites (n_seeds)
  ! LL: length of computational cube
  ! G: growing speed, scalar
  ! n_rnd_pos: number of random positions to check for crystallization time
  ! t_conv: to be returned, time of conversion at the positions (n_rnd_pos)
  ! vols: approximate volumes of the spherulites
  
  integer, intent(in) :: n_seeds, n_rnd_pos
  real, dimension(n_seeds, 3), intent(in) :: centers
  real, dimension(n_seeds), intent(in) :: t0s
  real, intent(in) :: LL, G
  real, dimension(n_rnd_pos), intent(out) :: t_conv
  real, dimension(n_seeds,2), intent(out) :: vols
  
  real, dimension(n_rnd_pos, 3) :: test_points
  real, dimension(n_seeds) :: t_conv_arr
  
  integer :: i, mloc

  real, dimension(n_seeds) :: dd_arr

  call init_random_seed()
  call random_number(test_points)
  test_points = test_points*LL
  vols(:,1) = 0.0
  
  do i = 1, n_rnd_pos
    call arr_distance(test_points(i,:), centers, n_seeds, dd_arr)
    t_conv_arr = dd_arr/G+t0s
    mloc = minloc(t_conv_arr,1)
    t_conv(i) = t_conv_arr(mloc)
    vols(mloc,1) = vols(mloc,1)+1.0


 end do

 vols(:,1) = vols(:,1)/n_rnd_pos*LL**3
 vols(:,2) = t0s
end subroutine sr_cdf2


!             _| |_               
!             \   /              
!              \ /               
!               Y  

subroutine sr_cdf2_aniso(n_seeds, centers, t0s, LL, T_0, Cr, tmax, &
                 n_rnd_pos, params, t_conv, vols)
  ! same as sr_cdf2 but for anisotherm case
  !  calculates the cdf of a single realization, given the positions and 
  ! starting time of the seeds. 
  ! The subroutine assumes that all the seeds are born. Thus please run 
  ! sp_seeds subroutine before this one, if calculating spontaneous crystal growth.
  !
  ! n_seeds: number of spherulite seeds to grow
  ! centers: coordinates of the centers of spherulites (n_seeds,3)
  ! t0s: birthtime of the spherulites (n_seeds)
  ! LL: length of computational cube
  ! T_0 - starting temperature of crystallization [K]
  ! Cr: cooling rate [K/s]
  ! params: fitted parameters of the Hoffmann-Lauritzen equation
  ! n_rnd_pos: number of random positions to check for crystallization time
  ! t_conv: to be returned, time of conversion at the positions (n_rnd_pos)
  ! vols: approximate volumes of the spherulites + their birthtimes
  
  
  
  
  integer, intent(in) :: n_seeds, n_rnd_pos
  real, dimension(n_seeds, 3), intent(in) :: centers
  real, dimension(n_seeds), intent(in) :: t0s
  real, intent(in) :: LL, T_0, Cr, tmax
  real, dimension(2), intent(in) :: params
  real, dimension(n_rnd_pos), intent(out) :: t_conv
  real, dimension(n_seeds,2), intent(out) :: vols
  
  real, dimension(n_rnd_pos, 3) :: test_points
  real, dimension(n_seeds) :: t_conv_arr
  
  integer :: i, mloc, j

  real, dimension(n_seeds) :: dd_arr, t1_arr
  real :: t_temp, ddmax
  real, dimension(3) :: g00

  call init_random_seed()
  call random_number(test_points)
  test_points = test_points*LL
  vols(:,1) = 0.0
  
  call HL_calc_ido(T_0, 0.0, Cr, params, g00(1))
  call G_elso(T_0, Cr, params, g00(2))
  call G_masodik(T_0, Cr, params, g00(3))
  call dist_calc(tmax, 0.0, T_0, Cr, params, ddmax)
  
  do i = 1, n_rnd_pos
    call arr_distance(test_points(i,:), centers, n_seeds, dd_arr)
    do j = 1, n_seeds
        if (dd_arr(j)<ddmax) then !megnézi, hogy a távolság nagyobb-e egy elméleti maximumnál, ha igen, akkor nem számolja ki az időt
        call time_calc(dd_arr(j), t0s(j), T_0, Cr, params, g00, t_temp) 
            t1_arr(j) = t_temp
        else
            t1_arr(j)=tmax+1
        endif
    end do
    t_conv_arr = t1_arr + t0s
    
    
    mloc = minloc(t_conv_arr,1)
    t_conv(i) = t_conv_arr(mloc)
    vols(mloc,1) = vols(mloc,1)+1.0


 end do

 vols(:,1) = vols(:,1)/n_rnd_pos*LL**3
 vols(:,2) = t0s
end subroutine sr_cdf2_aniso

!            / \     
!           /   \
!          /_   _\  
!            | | 




function distance(c1,c2)
  real :: distance
  real, dimension(3) :: c1, c2
  
  distance = sqrt(sum((c1-c2)**2))
  return
end function distance


subroutine arr_distance(sc, arr, arr_length, dd)

  integer, intent(in) :: arr_length
  real, dimension(arr_length, 3), intent(in) :: arr
  real, dimension(3), intent(in) :: sc
  real, dimension(arr_length), intent(out) :: dd
  
  dd = sqrt(sum((arr-spread(sc,1,arr_length))**2, 2))
  
end subroutine arr_distance


subroutine init_random_seed()

  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
end subroutine init_random_seed


subroutine MSort(N, A, A_out)
!a sorting algorithm, sorts A that is length N 
!returns A_out, sorthed array
        integer, intent(in) :: N 
        real, dimension(N), intent(in) :: A
        real, dimension(N), intent(out) :: A_out
        real :: work((size(A) + 1) / 2)
        A_out=A
        call MergeSort(A_out, work)
      contains
 
      subroutine merge(A, B, C)
! The targe attribute is necessary, because A .or. B might overlap with C.
        real, target, intent(in) :: A(:), B(:)
        real, target, intent(inout) :: C(:)
        integer :: i, j, k
        if (size(A) + size(B) > size(C)) then
            stop
        end if
        i = 1; j = 1
        do k = 1, size(C)
          if (i <= size(A) .and. j <= size(B)) then
            if (A(i) <= B(j)) then
              C(k) = A(i)
              i = i + 1
            else
              C(k) = B(j)
              j = j + 1
            end if
          else if (i <= size(A)) then
            C(k) = A(i)
            i = i + 1
          else if (j <= size(B)) then
            C(k) = B(j)
            j = j + 1
          end if
        end do
      end subroutine merge
 
      subroutine swap(x, y)
        real, intent(inout) :: x, y
        real :: tmp
        tmp = x; x = y; y = tmp
      end subroutine
 
      recursive subroutine MergeSort(A, work)
        real, intent(inout) :: A(:)
        real, intent(inout) :: work(:)
        integer :: half
        half = (size(A) + 1) / 2
        if (size(A) < 2) then
          continue
        else if (size(A) == 2) then
          if (A(1) > A(2)) then
            call swap(A(1), A(2))
          end if
        else
          call MergeSort(A( : half), work)
          call MergeSort(A(half + 1 :), work)
          if (A(half) > A(half + 1)) then
            work(1 : half) = A(1 : half)
            call merge(work(1 : half), A(half + 1:), A)
          endif
        end if
      end subroutine MergeSort
end subroutine MSort


!             _| |_               
!             \   /              
!              \ /               
!               Y       

subroutine cubic(g0, g1, g2, dd, real_root)
!cubic equation solver, returns one real root
    real, intent(in) :: g0, g1, g2, dd
    real, intent(out) :: real_root
    real :: dis, a, b, c, d
    real :: s, t, q, r, temp
    
    
    
    a = g2/6
    b = g1/2
    c = g0
    d = -dd
    dis = (18*a*b*c*d) - (4*b**3*d) + (b**2*c**2) - (4*a*c**3) - (27*a**2*d**2)
    q       = (3.0*a*c - b**2)/(9.0*a**2) 
    r       = (9.0*a*b*c - 27.0*d*a**2 - 2.0*b**3)/(54.0*a**3)

    temp    = r + sqrt(q**3 + r**2)
    s       = sign(1.0, temp) * abs(temp)**(1.0/3.0)

    temp    =  r - sqrt(q**3 + r**2)
    t       = sign(1.0, temp) * abs(temp)**(1.0/3.0)


    if (dis<0) then
        real_root = s + t - b/(3.0*a)
    end if


end subroutine cubic


subroutine HL_calc_ido(T0, ido, Cr, params, G)
    real, intent(in) :: ido, T0, Cr
    real, dimension(2), intent(in) :: params
    real, intent(out) :: G
    real :: G0, U, Tg, T0m, KG, Tref
    real :: T1, T2, T3, A1, A2, A3, A4, T
    !kiszámolja G értékét egy adott időpillanatban ha ismert a hűtési sebesség és a kezdeti hőmérséklet
    !calculates the growth speed at a certain time from the Hoffmann-Lauritzen equation
    !T0 : starting temperature of crystallization [K]
      ! Cr :cooling rate [K/s]
      ! params: fitted parameters of the Hoffmann-Lauritzen equation
      ! ido: the time at which we want to know G
      !returns G: growth speed [m/s]

    G0=params(1)
    U=755 
    Tg=263.15 
    Tref=30.0 
    T0m=481.15 
    KG=params(2) 
    
    T1 = -T0+Tg-Tref
    T2 = T0m+T0
    T3 = T0m-T0
    A1 = (KG*T0m**2)/2
    A2 = T3*T0**2
    A3 = Cr*(T0**2-2*T3*T0)
    A4 = Cr**2*(T3-2*T0)
    
    T = (T0-Cr*ido)

    G=G0*exp(-U/(T-Tg+Tref))*exp(KG*(T0m**2)*(T0m+T)/(2*(T**2)*(T0m-T)))


end subroutine HL_calc_ido



subroutine HL_calc(T, params, G)
    real, intent(in) :: T
    real, dimension(2), intent(in) :: params
    real, intent(out) :: G
    real :: G0, U, Tg, T0m, KG, Tref
    !calculates G at a given temperature from the Hoffmann-Lauritzen equation
    !params: fitted parameters of the Hoffmann-Lauritzen equation
    !T: temeprature [K]
    !returns G: growth speed [m/s]

    G0=params(1)
    U=755 
    Tg=263.15 
    Tref=30.0 
    T0m=481.15 
    KG=params(2) 
    
    G=G0*exp(-U/(T-Tg+Tref))*exp(KG*(T0m**2)*(T0m+T)/(2*(T**2)*(T0m-T)))

end subroutine HL_calc


subroutine G_elso(T0, Cr, params, eredmeny)
    real, intent(in) :: T0, Cr
    real, dimension(2), intent(in) :: params
    real, intent(out) :: firstderivate
    real :: G0, U, Tg, T0m, KG, Tref
    real :: T1, T2, T3, A1, A2, A3, A4
    !calculates first order derivate of G
    !T0 : starting temperature of crystallization [K]
      ! Cr :cooling rate [K/s]
      ! params: fitted parameters of the Hoffmann-Lauritzen equation
      ! returns: firstderivate

    G0=params(1) 
    U=755 
    Tg=263.15 
    Tref=30.0 
    T0m=481.15 
    KG=params(2) 
    
    T1 = -T0+Tg-Tref
    T2 = T0m+T0
    T3 = T0m-T0
    A1 = (KG*T0m**2)/2
    A2 = T3*T0**2
    A3 = Cr*(T0**2-2*T3*T0)
    A4 = Cr**2*(T3-2*T0)
    
    firstderivate = G0*(-Cr*A1/A2 - A1*T2*A3/A2**2 - U*Cr/T1**2)*exp(A1*T2/A2 + U/T1)

end subroutine G_elso


subroutine G_masodik(T0, Cr, params, eredmeny)
    real, intent(in) :: T0, Cr
    real, dimension(2), intent(in) :: params
    real, intent(out) :: secondderivate
    real :: G0, U, Tg, T0m, KG, Tref
    real :: T1, T2, T3, A1, A2, A3, A4
    real :: eredmeny1, eredmeny2
        !calculates fsecond order derivate of G
    !T0 : starting temperature of crystallization [K]
      ! Cr :cooling rate [K/s]
      ! params: fitted parameters of the Hoffmann-Lauritzen equation
      ! returns: secondderivate

    G0=params(1) 
    U=755 
    Tg=263.15 
    Tref=30.0 
    T0m=481.15 
    KG=params(2) 
    
    T1 = -T0+Tg-Tref

    T2 = T0m+T0
    T3 = T0m-T0
    A1 = (KG*T0m**2)/2
    A2 = T3*T0**2

    A3 = Cr*(T0**2-2*T3*T0)

    A4 = Cr**2*(T3-2*T0)

    
    eredmeny1 = ((-Cr*A1/A2 - A1*T2*A3/A2**2 - U*Cr/T1**2)**2)
    eredmeny2 = (2*Cr*A1*A3/A2**2 - A1*T2*2*A4/A2**2 + A1*T2*A3*2/A2**3 + 2*U*Cr**2/T1**3)
    secondderivate = G0*(eredmeny1 + eredmeny2)*exp(A1*T2/A2 + U/T1)

end subroutine G_masodik

subroutine G_harmadik(T0, Cr, params, eredmeny)
    real, intent(in) :: T0, Cr
    real, dimension(2), intent(in) :: params
    real, intent(out) :: eredmeny
    real :: G0, U, Tg, T0m, KG, Tref
    real :: T1, T2, T3, A1, A2, A3, A4
    real :: e1, e2, e3, e4, e5
    !G harmadik deriváltja

    G0=params(1) 
    U=755 
    Tg=263.15 
    Tref=30.0 
    T0m=481.15 
    KG=params(2)  
    
    T1 = -T0+Tg-Tref
    T2 = T0m+T0
    T3 = T0m-T0
    A1 = (KG*T0m**2)/2
    A2 = T3*T0**2
    A3 = Cr*(T0**2-2*T3*T0)
    A4 = Cr**2*(T3-2*T0)
    
    e1 = (-Cr*A1/A2 - A1*T2*A3/A2**2 - U*Cr/T1**2)**3
    e2 = 2*Cr*A1*A3/A2**2 - A1*T2*2*A4/A2**2 + A1*T2*A3*2/A2**3 + 2*U*Cr**2/T1**3
    e3 = -Cr*A1/A2 - A1*T2*A3/A2**2 - U*Cr/T1**2
    e4 = ((Cr*A1*A4-Cr**3*A1*T2)/A2**2 - (Cr*A1*A3**2 + A1*T2*2*A4*A3)/A2**3)
    e5 = A1*T2*A3**3/A2**4 - U*Cr**3/T1**4
    eredmeny = G0*(e1 + 3*e2*e3+6*(e4+e5))*exp(A1*T2/A2 + U/T1)

end subroutine G_harmadik

subroutine time_calc(d, t0, T_0, Cr, params, g00, t1) 
    real, intent(in) :: d, t0, T_0, Cr
    real, dimension(3), intent(in) :: g00 !T0=0 esetre
    real, dimension(2), intent(in) :: params
    real, intent(out) :: t1 
    real :: G0, G1, G2, dd
    !calculates how many seconds it takes for the crystal to grow d distance
    !d: distance [um/s]
    !T_0 : starting temperature of crystallization [K]
      ! Cr :cooling rate [K/s]
      ! params: fitted parameters of the Hoffmann-Lauritzen equation

    !t0: time of starting to grow [s]
    !g00: if t0=0, then these are the derivates [array of 3] (only needed for optimization)
    !returns: t1-the time it takes to grow [s]

    dd = d/(10**6)

    if (T_0==0) then
        G0 = g00(1)
        G1 = g00(2)
        G2 = g00(3)
    else
        call HL_calc_ido(T_0, t0, Cr, params, G0)
        call G_elso(T_0, Cr, params, G1)
   ! write(*,*)'első derivált :',G1
        call G_masodik(T_0, Cr, params, G2)
   ! write(*,*)'második derivált :',G2

    end if
    call cubic(G0, G1, G2, dd, t1)
end subroutine time_calc


subroutine dist_calc(t1, t0, T_0, Cr, params, d)
    real, intent(in) :: t1, t0, T_0, Cr
    real, dimension(2), intent(in) :: params
    real, intent(out) :: d
    real :: G0, G1, G2
    !calculates how much the crystal grows in t1 seconds starting from t0 time
    !returns d: distance [um/s]
    !T_0 : starting temperature of crystallization [K]
      ! Cr :cooling rate [K/s]
      ! params: fitted parameters of the Hoffmann-Lauritzen equation

    !t0: time of starting to grow [s]
    !t1-the time it takes to grow [s]
       

    call HL_calc_ido(T_0, t0, Cr, params, G0)
    call G_elso(T_0, Cr, params, G1)
    call G_masodik(T_0, Cr, params, G2)
    
    d = (G0*t1 + G1*t1**2/2 + G2*t1**3/6)*10**6 !m i k r o m é t e r

end subroutine dist_calc