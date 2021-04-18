      PROGRAM bene

* ---- variabili di sistema

      implicit none

      CHARACTER*16 nome
      CHARACTER*70 mess
      INTEGER*4 len, n, i1, j, k, l, ls
      REAL*8 p0d_v(17), t0_v(1)
      REAL*8 p0d, t0
      REAL*8 k_b, msuk, gamma, r_l, h_0, pas, acc
      REAL*8 M, t_l0, t_l, u_l, n_l, fi
      REAL*8 dat_T(100000), dat_O(100000), c(100000,3)
      REAL*8 cd(3), cd_old(3), x(3), y(3)
      REAL*8 r, n_r, u_r, t_r(2), tp_r(2), t_r1(2), t_rn(2)
      REAL*8 sr, h, rapp
      LOGICAL flag

      COMMON / integ / c, dat_T, len
      COMMON / par / u_r, n_r

      parameter (k_b = 1.380662d-23, msuk = 4.814053d-4)
      parameter (r_l = 2.5, gamma = 5.0 / 3.0)
      parameter (h_0 = 1.d-3, pas = 1.001, acc = 1.d-2)

* ---- definisco i valori di p0d e t0 -> defining p0d e t0

      DATA  p0d_v /1,2,3,5,8,10,15,20,30,50,81,100,120,130,150,180,200/

      DATA t0_v /118.5/

* ---- leggo i valori di omega(T) -> reading values of omega(T)

      open (8, file = 'omega_LJ.dat', err = 15, status = 'old')
      len = 0
      do n = 1, 100001
          read (8, fmt = *, err = 20, end = 200) dat_T(n), dat_O(n)
          len = len +1
      enddo

      mess = 'ATTENZIONE !!! Non tutti i valori di omega(T) sono ' //
     1        'stati caricati. Error reading'
      goto 210

 200  mess = 'O.K. valori di omega(T) letti correttamente !!'

 210  write(*,*) mess
      close(8)

* ---- interpola omega al variare di T con un polinomio del secondo ordine-> second order interpolation

      write(*,*) 'Interpolo omega(T) .....'
      do i1 = 2, len-1
          do n = 1, 3
              cd_old(n) = cd(n)
              x(n) = dat_T(i1-2+n)
              y(n) = dat_O(i1-2+n)
          enddo
          call terzo(x, y, cd)
          if (i1 .eq. 2) then
              do n = 1, 3
                  c(i1-1, n) = cd(n)
              enddo
          else
              do n = 1, 3
                  c(i1-1, n) = (cd(n) + cd_old(n)) / 2
                  if (i1 .eq. len-1) c(i1, n) = cd(n)
              enddo
          endif
      enddo
      write(*,*) 'Fatto !!! Done!!!'

      do l = 1, 1 !temperature
      do k = 1, 17 !pressure (really p0d but with 10 micron nozzle=10^-3 cm and pressure in bar -> unit is mbar*cm)

* ---- costruisco il nome del file di uscita -> constraction output file name (nome)

          t0 = t0_v(l)
          p0d = p0d_v(k)
          nome(1:1) = 't'
          write(unit = nome(2:4), fmt = '(i3)') nint(t0)
          nome(5:5) = 'p'
          if ((p0d .ge. 1) .and. (p0d .lt. 10)) then
              write(unit = nome(6:6), fmt = '(i1)') nint(p0d)
              ls = 6
          elseif ((p0d .ge. 10) .and. (p0d .lt. 100)) then
              write(unit = nome(6:7), fmt = '(i2)') nint(p0d_v(k))
              ls = 7
          elseif ((p0d .ge. 100) .and. (p0d .lt. 1000)) then
              write(unit = nome(6:8), fmt = '(i3)') nint(p0d_v(k))
              ls = 8
          else
              write(*,*) 'Errore nella scrittura del nome del file ' //
     1                   'di output !!'
          endif
          nome(ls+1 : ls+2) = 'LJ'
          nome(ls+3 : ls+6) = '.dat'

          write(*,*) 'p0d = ', p0d, ' t0 = ', t0
          write(*,*) 'file di uscita -> ', nome(1:ls+8)
          write(*,*)'Calcolo le condizioni iniziali .....'

* ---- condizioni iniziali (r_l = 2.5) in approssimazione sferica -> initial conditions in spherical approx.

          M = ( r_l**(gamma-1.d0) )*(3.232 - .7563/r_l +
     1        .3937/(r_l**2) - .0729/(r_l**3))
          t_l0 = 1 / (1 + (gamma-1.d0)/2 * M**2)
          t_l = t0 * t_l0
          u_l = M * sqrt(gamma * t_l / msuk)
          n_l = p0d / (k_b*t0) * t_l0**(1.d0 / (gamma - 1.d0))

          open (13, file = nome(1 : ls+6), err = 25)
          write (13, 100, err = 30) r_l, n_l, u_l, t_l, t_l
          close (13)

          write(*,*) 'Fatto !!! Done!!!'

* ---- risoluzione sistema -> Solving the differential equations

          write(*,*) 'Risolvo il sistema di equazioni differenziali ...'

          fi = n_l * u_l * r_l**2
          r = r_l
          u_r = u_l
          n_r = n_l
          t_r(1) = t_l
          t_r(2) = t_l
          j = 0

          flag = .true.
          do while (flag)

              h = h_0 * pas**j

              call derivate(r, t_r, tp_r)
              t_rn = t_r + h * tp_r / 6.

              r = r + h/2.
              do i1 = 1, 2
                  t_r1 = t_r + h/2. * tp_r
                  call derivate(r, t_r1, tp_r)
                  t_rn = t_rn + h * tp_r / 3.
              enddo

              r = r + h/2.
              t_r1 = t_r + h * tp_r
              call derivate(r, t_r, tp_r)
              t_r = t_rn + h * tp_r / 6.

              u_r = sqrt((5.d0*t0 - 3.d0*t_r(2) - 2.d0*t_r(1)) / msuk)
              n_r = fi / (u_r * r**2)

              open (13, file = nome(1 : ls+6), access = 'append',
     1              err = 25)
              write (13, 100, err = 30) r, n_r, u_r, t_r
              close (13)

              rapp = t_r(1) / t_r(2)
              if (rapp .lt. acc) then
                  flag = .false.
                  sr = sqrt(msuk * u_r**2. / (2. * t_r(2)))
                  write(*,*) 'Speed Ratio = ', sr
                  open (14, file = 'sp_raLJ_He.dat', access = 'append',
     1                  err = 35)
                  write (14, 101, err = 40) p0d, t0, sr
                  close (14)
              endif
              j = j + 1
          enddo
      enddo
      enddo

      write(*,*) 'O.K.!!! programma terminato correttamente'

      goto 1000

 100  format (f12.4, 1x, e10.4, f10.4, 1x, e10.4, 1x, e10.4)
 101  format (3f20.5)


* ---- gestione errori -> error handling

 10   mess = 'Errore nella apertura del file contenente p0d e T0 !!'
      goto 50
 15   mess = 'Errore nella apertura del file contenente omega(T) !!'
      goto 50
 20   mess = 'Errore nella lettura del file contenente omega(T) !!'
      goto 50
 25   mess = 'Errore nella apertura del file 1 di output !!'
      goto 50
 30   mess = 'Errore nella scrittura del file 1 di output !!'
      goto 50
 35   mess = 'Errore nella apertura del file 2 di output !!'
      goto 50
 40   mess = 'Errore nella scrittura del file 2 di output !!'
      goto 50


 50   write(*,*) char(7)
      write(*,*) mess
      stop


 1000 end




* ---- calcolo le derivate prime di T_par e T_per in r -> derivative calculations

      SUBROUTINE derivate(r, t_r, tp_r)

      implicit none

      REAL*8 r, t_r(2), tp_r(2), ris, t1
      REAL*8 u_r, n_r, msuk

      COMMON / par / u_r, n_r
      PARAMETER(msuk = 4.814053d-4)

      call integrale(t_r, ris)

      ris = 2 * (n_r / u_r) * ris / (t_r(1) * sqrt(t_r(2)))

      t1 = u_r**2. * msuk
      tp_r(1) = -2.d0 * t_r(1) / r + ris
      tp_r(2) = (-2.d0 * ris * t1 + 2.d0 * t_r(2) *
     1          (-2.d0 * t_r(1) / r + ris)) / (t1 - 3.d0 * t_r(2))

      return
      end




* ---- calcolo l'integrale I -> integral calculation

      SUBROUTINE integrale(t_r, ris)

      implicit none

      REAL*8 c(100000, 3), dat_T(100000)
      REAL*8 t_r(2), ris
      REAL*8 chi, dchi, a, t_eff, omega, fun
      REAL*8 A_O
      INTEGER*4 i, len
      LOGICAL flag1, flag2

      COMMON / integ / c, dat_T, len

      parameter(A_O = 2.48004d-20)

      ris = 0.d0
      dchi = 1.d-2
      chi = dchi
      flag1 = .true.

      do while (flag1)

          a = (t_r(2) - t_r(1)) / t_r(2)
          t_eff = t_r(1) / (1 - a * chi**2.)
          flag2 = .true.
          i = 0

          do while (flag2)
              i = i + 1
              if (dat_T(i) .ge. t_eff) flag2 = .false.
          enddo

          if (i .eq. 1) then
              omega = 0
          elseif ((i .gt. 1) .and. (i .le. len)) then
              omega = A_O * ( c(i-1, 3) * t_eff**2. + c(i-1, 2)
     1                * t_eff + c(i-1, 1))
          else
              write(*,*) 'ERRORE NELLA ROUTINE DI INTEGRAZIONE !!'
              stop
          endif

          fun = t_eff**(5./2.) * omega * (3.* chi**2. - 1.)
          ris = ris + fun * dchi
          chi = chi + dchi
          if (chi .gt. 1) flag1 = .false.

      enddo

      return
      end




* ---- calcola i coefficienti del polinomio di secondo grado che interpola
* ---- omega(T) -> second order interpolation

      SUBROUTINE terzo(x, y, c)

      implicit none

      REAL*8 x(3), y(3), c(3)
      REAL*8 x21, x31, xq12, xq13, y21, y31, p, q

      x21 = x(2) - x(1)
      xq12 = x(1)**2 - x(2)**2
      x31 = x(3) - x(1)
      xq13 = x(1)**2 - x(3)**2
      y21 = y(2) - y(1)
      y31 = y(3) - y(1)
      p = y31 / x31 - y21 / x21
      q = xq12 / x21 - xq13 / x31

      c(3) = p / q
      c(2) = (c(3) * xq13 + y31) / x31
      c(1) = y(1) - c(3) * x(1)**2 - c(2) * x(1)

      return
      end
