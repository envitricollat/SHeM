      PROGRAM bene

* ---- variabili di sistema

      implicit none

      CHARACTER*16 nome
      CHARACTER*70 mess
      INTEGER*4 len, n, i1,j, k, l, ls, p0d_v(25), inis,n_press,n_temp
      REAL*8 p0d, t0
      REAL*8 k_b, msuk, gamma, r_l, h_0, pas, acc
      REAL*8 M, t_l0, t_l, u_l, n_l, fi
      REAL*8 dat_T(100000), dat_O(100000), c(100000,3), t0_v(2)
      REAL*8 cd(3), cd_old(3), x(3), y(3)
      REAL*8 r, n_r, u_r, t_r(2), tp_r(2), t_r1(2), t_rn(2)
      REAL*8 sr, h, rapp, T_Ent_E, rho_r, rho_real, dE_real
      LOGICAL flag

      COMMON / integ / c, dat_T, len
      COMMON / par / u_r, n_r

      parameter (k_b = 1.380662d-23, msuk = 4.814053d-4)
      parameter (r_l = 2.5, gamma = 5.0 / 3.0)
      parameter (h_0 = 1.d-3, pas = 1.001, acc = 1.d-2)

* ---- definisco i valori di p0d e t0

      DATA  p0d_v
     & /1,2,4,6,8,12,16,20,25,30,35,40,50,60,70,80,90,100,
     & 110,120,130,150,160,180,200/

C      DATA t0_v /300/
      DATA t0_v /131,132/
      n_press=25 !dimensione di p0d_v usato nel loop -> dimension of p0d_v used in the loop
      n_temp=2 !dimensione di t0_v usato nel loop -> dimension of t0_v used in the loop
* ---- leggo i valori di omega(T)

      open (8, file = 'omega_LJ.dat', err = 15, status = 'old')
      len = 0
      do n = 1, 100001
          read (8, fmt = *, err = 20, end = 200) dat_T(n), dat_O(n)
          len = len +1
      enddo

      mess = 'ATTENZIONE !!! Non tutti i valori di omega(T) sono ' //
     1        'stati caricati. '
      goto 210

 200  mess = 'O.K. valori di omega(T) letti correttamente !!'

 210  write(*,*) mess
      close(8)

* ---- interpola omega al variare di T con un polinomio del secondo ordine

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
      write(*,*) 'Fatto !!!'

      do l = 1, n_temp !temperature
      do k = 1, n_press !pressure

* ---- costruisco il nome del file di uscita

          t0 = t0_v(l)
          p0d = p0d_v(k)
!!!!!!!!!!!!!!!!correzione gas reale- vale solo per elio!!!!!!!! real gas correction valid only for He
          rho_r=rho_real(t0,p0d)
C         dE_real=2d0*p0d/rho_r/9.86896d-3/6.02214179d23/k_b !litri*atm in Joule/molecola from liter*bar to joule/molecule
          dE_real=2d0*p0d/rho_r/1d-2/6.02214179d23/k_b !litri*bar in Joule/molecola
          write(6,*)'Real gas correction: PV=',dE_real,' rho=',rho_r
          T_Ent_E=3.0d0*t0+dE_real !solo energia interna 3 kT+ PV ->internal energy 3 kT plus PV
!!!!!!!!!!!!!!!!

          nome(1:1) = 't'
        if(t0.lt.100)then
               write(unit = nome(2:3), fmt = '(i2)') nint(t0)
             inis=4
        else
               write(unit = nome(2:4), fmt = '(i3)') nint(t0)
             inis=5
        end if
          nome(inis:inis) = 'p'
          if ((p0d .ge. 1) .and. (p0d .lt. 10)) then
              write(unit = nome(inis+1:inis+1), fmt = '(i1)') nint(p0d)
              ls = inis+1
          elseif ((p0d .ge. 10) .and. (p0d .lt. 100)) then
              write(unit = nome(inis+1:inis+2), fmt = '(i2)') p0d_v(k)
              ls = inis+2
          elseif ((p0d .ge. 100) .and. (p0d .lt. 1000)) then
              write(unit = nome(inis+1:inis+3), fmt = '(i3)') p0d_v(k)
              ls = inis+3
          else
              write(*,*) 'Errore nella scrittura del nome del file ' //
     1                   'di output !!'
          endif
          nome(ls+1 : ls+2) = 'LJ'
          nome(ls+3 : ls+6) = '.dat'

          write(*,*) 'p0d = ', p0d, ' t0 = ', t0
          write(*,*) 'file di uscita -> ', nome(1:ls+8)
          write(*,*)'Calcolo le condizioni iniziali .....'

* ---- condizioni iniziali (r_l = 2.5) in approssimazione sferica

          M = ( r_l**(gamma-1.d0) )*(3.232 - .7563/r_l +
     1        .3937/(r_l**2) - .0729/(r_l**3))
          t_l0 = 1 / (1 + (gamma-1.d0)/2 * M**2)
          t_l = t0 * t_l0
          u_l = M * sqrt(gamma * t_l / msuk)
          n_l = p0d / (k_b*t0) * t_l0**(1.d0 / (gamma - 1.d0))

          open (13, file = nome(1 : ls+6), err = 25)
          write (13, 100, err = 30) r_l, n_l, u_l, t_l, t_l
          close (13)

          write(*,*) 'Fatto !!!'

* ---- risoluzione sistema

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

C              u_r = sqrt((5.d0*t0 - 3.d0*t_r(2) - 2.d0*t_r(1)) / msuk)
              u_r = sqrt((T_Ent_E - 3.d0*t_r(2) - 2.d0*t_r(1)) / msuk)
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
                  open (14, file = 'sp_raLJ.dat', access = 'append',
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


* ---- gestione errori

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




* ---- calcolo le derivate prime di T_par e T_per in r

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




* ---- calcolo l'integrale I

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
* ---- omega(T)

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
C------------

       real*8 function EOS(temp,rho)
C       real*8 function EOS_He(temp,dens)
!printf( "  Temperature   $temp  K ");
!printf( "  Density  $dens  kg/m3 ");  mol/L*amu=kg/m**3 so mol/L*4.0026
!pressure in bar
!The validity range of Hands and McCarty fits (used here)
!extends from 2 K up to 1500 K, for pressure up to 2 GPa
!Limited accuracy in the critical region (~ 5.19 K, ~ 69 kg/m3).
       implicit none
       real*8 rho, dens, temp, tempo1, tempo2, tempo3 !, pression
       real*8 tp2, tp3, tp4
       dens=rho*4.0026 !mol/L-> kg/m**3
       tp2=temp**2
       tp3=temp**3
       tp4=temp**4

       tempo2=2.367205158795324D-7*dens**11*(6.872567403738D-15/tp3
     1  - 4.595138561035D-15/tp2) + 1.477581743701978D-8*dens**13*
     2  (3.848665703556D-18/tp4 - 7.636186157005D-18/tp3
     3  -6.097223119177D-19/tp2)

        tempo3=0.00006075816692925613D0*dens**7 *
     1  (3.850153114958D-8/tp3 - 2.067693644675999D-8/tp2)+
     2  (3.792453641033501D-6)*dens**9*(-((1.888462892389D-12)/tp4)
     3  - 1.399040627D-11/tp2)+tempo2

        tempo1=(0.01559457081650663D0*dens**3
     1  * (0.009166109232806D0/tp3
     2  - 0.006885401367690001D0/tp2) + 0.000973394851465435D0*
     3  dens**5*(-(0.00003315398880031D0/tp4)- 6.544314242937D-6/tp2)
     4  +tempo3)

        tempo3 = 0.003896110232475551D0*dens**4*(2.6997269279D-6
     1  - 0.00003954146691114/temp - 5.093547838380999D-9*temp)

        tempo3=tempo3 + 0.01559457081650663D0*dens**3*
     1  (-0.00004708238429298D0 + 0.002410763742104D0/(tp2)
     2  + 0.001132915232587D0/temp +1.454229259623D-6*temp)

         tempo3=tempo3 + 0.06241882915014949D0*dens**2 *
     1  (-0.007139657549318D0 - 0.01589302471561998D0/tp2
     2  + 0.009728903861441D0/temp
     3  + 0.001260692007853D0*temp**0.5D0
     4  + 0.00004558980227431D0*temp)

       tempo3=tempo3 - (1.348439408105817D-18*dens**9)/tp2
     1  -(6.304713842604079D-15*dens**7)/temp
     2  + 0.002077227302253535D0*dens*temp

        EOS=10*(1.510671273545713D-12*dens**5 +
     1  tempo1/dexp(0.0002061897349793637D0*dens**2)
     2  + 0.00001517967494360068D0*dens**8*
     3  (3.298960057070999D-11/tp2 + 6.446881346447997D-13/temp)
     4  + 0.0002431906389510405D0*dens**6*(-(5.501158366750001D-8/tp2)
     5  + 1.050712335784999D-8/temp) + tempo3)


C *   McCarty and Arp (1990) Advances in cryogenic engineering.
C R. W. Fast. New York, Plenum Press. 35: 1465-1475.
C **  Echelle provisoire de température de 1976. Bureau
C international des poids et mesures, F-92310 Sevres, France
C *** Hands, B. A. and V. D. Arp (1981) A correlation of thermal
C conductivity data for helium. Cryogenics 21: 697-703.
       return
       end
C---------

       real*8 function rho_real(T,Pb)
       implicit none
       real*8 T,P,rho,dr,Pc,Pb  !Pb in bar
       integer*4 ind,indd
       real*8 EOS
       external EOS
C      P=Pb/1.01325 !from bar to atm
       P=Pb !from bar to bar
       rho=0.1d0*P
       dr=rho*1d-2
       ind=0
       indd=0
10     if(ind.eq.0) then
          Pc=EOS(T,rho)
          if(Pc.lt.P)then
            rho=rho+dr
            indd=0
          else
            if(indd.eq.0)then
              dr=0.25*dr
              indd=1
            endif
            rho=rho-dr
            if(dr.lt.1d-10)ind=1
          endif
          goto 10
       endif
       rho_real=rho
       return
       end
