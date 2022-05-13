!Data-driven and constrained optimization of semi-local exchange and
!non-local correlation functionals for materials and surface chemistry
!
!Kai Trepte and Johannes Voss
!ktrepte@slac.stanford.edu and vossj@slac.stanford.edu
!
!Journal reference: J. Comput. Chem. 43, 1104 (2022).
!https://doi.org/10.1002/jcc.26872
!
!
!
!Code for computing the VCML exchange energy and functional derivatives
!
!
!
!FORTRAN 90 subroutine for calculating the VCML exchange energy below
!
!For semi-local correlation, the VCML-rVV10 meta-GGA XC functional uses
!the same correlation as the MS0, MS1, and MS2 functionals, which is
!based on PBE correlation with the constant parameter
!beta=0.06672455060314922 replaced by a Wigner-Seitz radius rs dependent
!function [J. Sun, B. Xiao, and A. Ruzsinszky, J. Chem. Phys. 137,
!051101 (2012); J.P. Perdew, A. Ruzsinszky, G.I. Csonka, L.A. Constantin,
!J. Sun, Phys. Rev. Lett. 103, 026403 (2009).]
!beta(rs) = 0.06672455060314922 * (1.0 + 0.1*rs) / (1.0 + 0.1778*rs)
!Correspondingly modified routines for GGA correlation can be found in
!libxc with name 'GGA_C_REGTPSS'.
![S. Lehtola, C. Steigemann, M.J.T. Oliveira, and M.A.L. Marques;
!Software X 7, 1 (2018).]
!https://www.tddft.org/programs/libxc
!
!For non-local correlation, the VCML-rVV10 meta-GGA XC functional uses
!the rVV10 functional [O.A. Vydrov and T. Van Voorhis, J. Chem. Phys.
!133, 244103 (2010); R. Sabatini, T. Gorni, and S. de Gironcoli, Phys.
!Rev. B, 87, 041108 (2013)] with an optimized b-parameter of 15.35.


module vcmlconst
!VCML exchange enhancement expansion coefficients cij
double precision, parameter, dimension(64) :: vcml_cij = (/ &
&   1.1050362267560025d0, -0.1304673327239498d0, -0.25273044468938444d0, &
&   0.0020345583050872945d0, 0.009705556829333915d0, -0.0018727613481398786d0, &
&   0.005056319358478653d0, -0.0014994572626212954d0, 0.19526954394443446d0, &
&   0.12131628073942294d0, -0.013135604251829597d0, -0.016823429546012295d0, &
&   -0.0021100890252897446d0, -0.0016609256494831233d0, 0.0028206838819829017d0, &
&   0.00017309630990864668d0, -0.00068200282327089d0, 0.0012341314639045392d0, &
&   -0.000835331263170036d0, -7.823588139015819d-05, -0.0014878680171769923d0, &
&   0.005061925051098745d0, -0.007631605623646023d0, -0.01006770315965861d0, &
&   -0.00217177716567727d0, 0.0024977311122498513d0, -0.0008670535705479461d0, &
&   0.0027822064319562786d0, -0.0002571281595426713d0, -3.656012084198544d-05, &
&   -0.009195715678311926d0, 0.010726279571787276d0, -0.001432652476750007d0, &
&   0.0050995906979556666d0, 0.0003180493235941731d0, -0.004704436332280876d0, &
&   0.0009891355730978566d0, -0.0010249162124576494d0, 0.0008367073496483024d0, &
&   -0.00031389079758955066d0, -0.004500541251076788d0, 0.0016437722411542371d0, &
&   8.482767148525194d-05, -0.00019375881298946268d0, -7.261106354828029d-05, &
&   -0.0038541498256550073d0, -0.0031296536914037784d0, 0.0038758929812102785d0, &
&   0.00030574929164576756d0, 0.0005970286163074767d0, -0.0009048853909642742d0, &
&   -0.000689695394243961d0, 0.0001331797359718674d0, -0.007555456486598222d0, &
&   0.001864317026752979d0, -0.00019095139973664826d0, -0.002025317083565653d0, &
&   0.0023160016166370034d0, 0.00018939021743243079d0, 0.0004308565933608885d0, &
&   -1.792697304428732d-05, -0.0005194058669188706d0, -0.00018156466410673526d0, &
&   -0.00029476504977320184d0 /)

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
!vcml_exchange computes VCML exchange energy density (output: ex)
!and the derivatives d(ex*rho)/drho (output: drho),
!d(ex*rho)/dgrad (output: dgrad) and d(ex*rho)/dtau (output: dtau)
!
subroutine vcml_exchange(rho, grad, tau, ex, drho, dgrad, dtau)
    !
    !get VCML exchange enhancement expansion coefficients vcml_cij
    !and other useful constants from module 'vcmlconst' above
    use vcmlconst
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
    !VCML exchange enhancement
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
    !VCML exchange enhancement
    call legendre(Lalpha, dLalpha, hat_alpha)

    !compute VCML exchange enhancement
    Fx = 0d0
    dFs = 0d0
    dFalpha = 0d0
    k = 1
    do j = 1, 8
        do i = 1, 8
            Fx = Fx + vcml_cij(k) * Ls(i)*Lalpha(j)
            dFs = dFs + vcml_cij(k) * dLs(i)*Lalpha(j)
            dFalpha = dFalpha + vcml_cij(k) * Ls(i)*dLalpha(j)
            k = k + 1
        enddo
    enddo

    dFs = dFs * dhat_s
    dFalpha = dFalpha * dhat_alpha

    !compute exchange energy and derivatives wrt. rho, grad, and tau
    !Dirac (LDA) exchange energy density
    sx = diracx * rho13
    !VCML exchange energy density
    ex = sx * Fx
    !derivative of exchange energy wrt. rho
    drho = sx * (4d0/3d0*(Fx-s*dFs) + dFalpha*(tauw/tautf-5d0/3d0*alpha))
    !... wrt. grad
    dgrad = sx * (dFs*dsdgrad*rho - 0.25d0*dFalpha*grad/tautf)
    !... wrt. tau
    dtau = sx * dFalpha*rho/tautf
    !
end subroutine
