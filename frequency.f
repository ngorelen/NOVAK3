      subroutine frequency(gamma_m,gamma_k,omega_s,omega_t,beta,omega,ik
     & )
c
c     To calculate the frequencies (omega(2)) as a solutions of the
c     equation: 
c     \[
c     -i \omega (1+\omega_T^2/(\omega_s^2-\omega^2))=
c                                               =\gamma_m+\beta*\gamma_k,
c     \]
c     where $\gamma_k$ and $\omega $ are of a complex type.
c
c     The solutions are:
c     \[
c     \omega=+-\sqrt{(\omega_s^2+\omega_T^2-\gamma^2)+-
c     \sqrt{(\omega_s^2+\omega_T^2-\gamma^2)^2+4\gamma^2\omega_s^2}}.
c     \]
c     Only one pair of this solutions satisfies that equation.
c
c     On the Output are:
c     omega - two dimensional complex array of frequencies;
c     ik    - the integer number of solutions. Usually ik equals 2, but 
c              ik may be equal 1, if \omega_T=0.
c
      complex gamma_k,omega(2),gamma,gamma2,a
      real gamma_m,omega_s,omega_t,beta
      do i=1,2
         omega(i)=0.
      enddo
      gamma=gamma_m+beta*gamma_k
      if(omega_t.eq.0.) then
         ik=1
         omega(1)=(0.,1.)*gamma
         return
      endif
      ik=2
      omega_t2=omega_t**2
      omega_s2=omega_s**2
      gamma2=gamma**2
      a=omega_s2+omega_t2-gamma2
      omega(1)=csqrt((a+csqrt(a**2+4.*gamma2*omega_s2))/2.)
      omega(2)=csqrt((a-csqrt(a**2+4.*gamma2*omega_s2))/2.)
      do i=1,2
         a=(0.,-1.)*omega(i)*csqrt(1.+omega_t2/(omega_s2-omega(i)**2))
         if(real(a)*real(gamma).lt.0..or.aimag(a)*aimag(gamma).lt.0.)
     &        then
            omega(i)=-omega(i)
         endif
      enddo
      return
      end
c
c**************************************************************************
