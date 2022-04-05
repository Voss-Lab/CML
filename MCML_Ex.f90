!MCML: combining physical constraints with experimental data for
!a multi-purpose meta-generalized gradient approximation
!
!Kris Brown, Yasheng Maimaiti, Kai Trepte, Thomas Bligaard, and Johannes Voss
!vossj@slac.stanford.edu
!
!Journal reference: J. Comput. Chem. 42, 2004 (2021)
!https://doi.org/10.1002/jcc.26732
!
!
!
!Code for computing the MCML exchange energy and functional derivatives
!
!
!
!FORTRAN 90 subroutine for calculating the MCML exchange energy below
!
!For correlation, the MCML meta-GGA XC functional uses the same correlation
!as the MS0, MS1, and MS2 functionals, which is based on PBE correlation
!with the constant parameter beta=0.06672455060314922 replaced by a
!Wigner-Seitz radius rs dependent function
![J. Sun, B. Xiao, and A. Ruzsinszky, J. Chem. Phys. 137, 051101 (2012).
!J.P. Perdew, A. Ruzsinszky, G.I. Csonka, L.A. Constantin, J. Sun;
!Phys. Rev. Lett. 103, 026403 (2009).]
!beta(rs) = 0.06672455060314922 * (1.0 + 0.1*rs) / (1.0 + 0.1778*rs)
!
!Correspondingly modified routines for GGA correlation can be found in
!libxc with name 'GGA_C_REGTPSS'.
![S. Lehtola, C. Steigemann, M.J.T. Oliveira, and M.A.L. Marques;
!Software X 7, 1 (2018).]
!https://www.tddft.org/programs/libxc


module mcmlconst
!MCML exchange enhancement expansion coefficients cij
double precision, parameter, dimension(64) :: mcml_cij = (/ &
&   1.0678412675928217d0, -0.15872242282520396d0, -0.23727374470038587d0, &
&   0.0025587527432851254d0, 0.006748483298726394d0, 0.0011994362281622633d0, &
&   0.0015528466146458777d0, 0.000584892206996479d0, 0.20323990913830245d0, &
&   0.11793635648230213d0, -0.01437960658302686d0, -0.010305714294261083d0, &
&   0.016832150866862326d0, -0.00025773333827270803d0, 0.0023346167766491333d0, &
&   0.0003837976998664341d0, -0.0006952718706718514d0, 0.00179463855686441d0, &
&   -0.001153807045825489d0, -0.0007090296813211244d0, 0.00013702886354574704d0, &
&   0.006670848599065867d0, -0.005498112922165805d0, 0.0014213910238437613d0, &
&   -0.0025656924772691145d0, 0.002125332357775206d0, -0.0009641371299507833d0, &
&   0.0037127861713210433d0, 0.0012824718527707636d0, 0.0002262886186270548d0, &
&   -0.006510071882485726d0, 0.01243327883803539d0, -0.002170152177993684d0, &
&   0.002915285520983635d0, -0.0018638828810102481d0, -0.0024949505505474645d0, &
&   0.00041878279077109046d0, -0.0010099812635462266d0, 0.0004230264400260503d0, &
&   0.0024575259185362595d0, -0.0027233877043555677d0, 0.002007295399058147d0, &
&   -0.0011896683049514127d0, 0.0001672905908063297d0, -0.0002721968500889238d0, &
&   -0.0005869916483960576d0, 0.0011364858250944854d0, 0.0015224741795989722d0, &
&   -0.0019776072156133598d0, 0.001491587478361034d0, -0.0012883061272796169d0, &
&   -0.0006058496834176058d0, 0.0002776060240069905d0, -0.0016226213909532258d0, &
&   0.0004260858412001439d0, -0.00036825194324629357d0, -0.002516160322803815d0, &
&   0.0019401647142238965d0, -0.0011756144767584228d0, 0.00043124117592430525d0, &
&   3.212943141118693d-6, -0.0002202759704065197d0, 0.00038071585953508916d0, &
&   -0.0003695503801501715d0 /)

!eta = kappa/muGE = 0.804/(10/81)
double precision, parameter :: eta = 6.5124d0

!sconst = (3/pi)**(2/3) / 6
double precision, parameter :: sconst = 0.16162045967399548133166d0

!tfconst = 3/10 * (3*pi**2)**(2/3)
double precision, parameter :: tfconst = 2.87123400018819181594250d0

!diracx = -3/4 * (3/pi)**(1/3)
double precision, parameter :: diracx = -0.73855876638202240588423d0
!
end module


!calculate Legendre polynomials L up to order 7 and first derivatives dL at x
subroutine legendre(L, dL, x)
    implicit none
    double precision L(*), dL(*), x
    !local variale
    integer i

    L(1) = 1d0
    L(2) = x
    dL(1) = 0d0
    dL(2) = 1d0
    !use recursion to compute Legendre polynomials order 2 to 7
    do i = 2, 7
        L(i+1) = 2d0*x*L(i) - L(i-1) - (x*L(i) - L(i-1))/i
    enddo
    !same for derivatives
    do i = 2, 7
        dL(i+1) = L(i)*i + dL(i)*x
    enddo
end subroutine


!
!Based on input density rho, length of density gradient grad, kinetic
!energy density tau (all at point r),
!mcml_exchange computes MCML exchange energy density (output: ex)
!and the derivatives d(ex*rho)/drho (output: drho),
!d(ex*rho)/dgrad (output: dgrad) and d(ex*rho)/dtau (output: dtau)
!
subroutine mcml_exchange(rho, grad, tau, ex, drho, dgrad, dtau)
    !
    !get MCML exchange enhancement expansion coefficients mcml_cij
    !and other useful constants from module 'mcmlconst' above
    use mcmlconst
    implicit none
    double precision rho, grad, tau, ex, drho, dgrad, dtau
    !local variables
    double precision Ls(8), dLs(8), Lalpha(8), dLalpha(8)
    double precision rho13, rho43, dsdgrad, s, s2, etapluss2
    double precision tauw, tautf, alpha, alpha2, alpha3
    double precision oneminusalpha2, oneminusalpha2_2, idenom
    double precision hat_s, dhat_s, hat_alpha, dhat_alpha
    double precision Fx, dFs, dFalpha, sx
    integer i, j, k

    rho13 = rho**(1d0/3d0)
    rho43 = rho*rho13
    dsdgrad = sconst / rho43
    !reduced density gradient
    s = dsdgrad * grad
    s2 = s*s
    etapluss2 = eta + s2
    hat_s = 2d0*s2 / etapluss2 - 1d0
    dhat_s = 4d0*eta*s / (etapluss2*etapluss2)
    !evaluate Legendre polynomials for gradient dependence of
    !MCML exchange enhancement
    call legendre(Ls, dLs, hat_s)

    !tau_Weizsaecker
    tauw = 0.125d0 * grad*grad / rho
    !tau_UEG
    tautf = tfconst * rho43*rho13
    alpha = (tau - tauw) / tautf
    if(alpha<1d6) then
        alpha2 = alpha*alpha
        alpha3 = alpha*alpha2
        oneminusalpha2 = 1d0 - alpha2
        oneminusalpha2_2 = oneminusalpha2*oneminusalpha2
        idenom = 1d0 / (1d0 + alpha3 + 4d0*alpha3*alpha3)
        hat_alpha = oneminusalpha2*oneminusalpha2_2 * idenom
        dhat_alpha = -3d0*oneminusalpha2_2*alpha * idenom*idenom * &
&                    (2d0 + alpha + alpha3 + 8d0*alpha2*alpha2)
    else
        hat_alpha = -0.25d0
        dhat_alpha = 0d0
    endif
    !evaluate Legendre polynomials for kinetic energy density dependence of
    !MCML exchange enhancement
    call legendre(Lalpha, dLalpha, hat_alpha)

    !compute MCML exchange enhancement
    Fx = 0d0
    dFs = 0d0
    dFalpha = 0d0
    k = 1
    do j = 1, 8
        do i = 1, 8
            Fx = Fx + mcml_cij(k) * Ls(i)*Lalpha(j)
            dFs = dFs + mcml_cij(k) * dLs(i)*Lalpha(j)
            dFalpha = dFalpha + mcml_cij(k) * Ls(i)*dLalpha(j)
            k = k + 1
        enddo
    enddo

    dFs = dFs * dhat_s
    dFalpha = dFalpha * dhat_alpha

    !compute exchange energy and derivatives wrt. rho, grad, and tau
    !Dirac (LDA) exchange energy density
    sx = diracx * rho13
    !MCML exchange energy density
    ex = sx * Fx
    !derivative of exchange energy wrt. rho
    drho = sx * (4d0/3d0*(Fx-s*dFs) + dFalpha*(tauw/tautf-5d0/3d0*alpha))
    !... wrt. grad
    dgrad = sx * (dFs*dsdgrad*rho - 0.25d0*dFalpha*grad/tautf)
    !... wrt. tau
    dtau = sx * dFalpha*rho/tautf
    !
end subroutine
