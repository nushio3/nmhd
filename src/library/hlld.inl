/*
  Argument order is as follows
  1. common parameter,
  2. initial condition, in pairs
  3. output parameters
 */
__host__ __device__
void hlld_flux(Real Bx, 
	       
	       Real dens_L, Real dens_R, 
	       Real velx_L, Real velx_R, 
	       Real vely_L, Real vely_R, 
	       Real velz_L, Real velz_R, 
	       Real pres_L, Real pres_R, 
	       Real By_L,   Real By_R,   
	       Real Bz_L,   Real Bz_R,   

	       Real &Fdens, 
	       Real &Fmomx, Real &Fmomy, Real &Fmomz,
	       Real &Fetot, 
	       Real &F_By,  Real &F_Bz) {

  const Real TINY = 1.0e-10;

  Real signBx = 0;
  if (absR(Bx) > 0) signBx = Bx/absR(Bx);
  
  Real momx_L = dens_L*velx_L;
  Real momy_L = dens_L*vely_L;
  Real momz_L = dens_L*velz_L;
  
  Real momx_R = dens_R*velx_R;
  Real momy_R = dens_R*vely_R;
  Real momz_R = dens_R*velz_R;
  
  Real B2_L   = sq(Bx)     + sq(By_L)   + sq(Bz_L);
  Real v2_L   = sq(velx_L) + sq(vely_L) + sq(velz_L);
  Real etot_L = pres_L/(kGamma - 1) + 0.5*dens_L*v2_L + 0.5*B2_L;

  Real B2_R   = sq(Bx)     + sq(By_R)   + sq(Bz_R);
  Real v2_R   = sq(velx_R) + sq(vely_R) + sq(velz_R);
  Real etot_R = pres_R/(kGamma - 1) + 0.5*dens_R*v2_R + 0.5*B2_R;

  Real gpl  = kGamma * pres_L;
  Real gpr  = kGamma * pres_R;
  Real gpbl = gpl + B2_L;
  Real gpbr = gpr + B2_R;

  const float fsig = sqrtR(2.0); // slightly increase S_L & S_R to avoid problems
  
#if 1

  // use magnetosonic speed for left and right singal speeds

  Real cfl = sqrtR((gpbl + sqrtR( sq(gpbl) - 4*gpl*sq(Bx) ))/(2.0*dens_L));
  Real cfr = sqrtR((gpbr + sqrtR( sq(gpbr) - 4*gpr*sq(Bx) ))/(2.0*dens_R));
  Real cfmax = fsig*fmax(cfl,cfr);

  Real S_R, S_L;

  S_L = fmin(velx_L, velx_R) - cfmax;
  S_R = fmax(velx_L, velx_R) + cfmax;

#else

  // Use Roe-average
  
  Real dsl = sqrtR(dens_L);
  Real dsr = sqrtR(dens_R);
  Real ids = 1.0/(dsl + dsr);
  Real droe = dsl * dsr;
  
  Real uxroe = (dsl*velx_L + dsr*velx_R)*ids;
  Real uyroe = (dsl*vely_L + dsr*vely_R)*ids;
  Real uzroe = (dsl*velz_L + dsr*velz_R)*ids;
  
  Real byroe = (dsl*By_L + dsr*By_R)*ids;
  Real bzroe = (dsl*Bz_L + dsr*Bz_R)*ids;
  
  Real x = 0.5 * (sq(By_L - By_R) + sq(Bz_L - Bz_R))/sq(dsl + dsl);
  Real y = 0.5 * (dens_L + dens_R)/droe;
    
  Real pbl = 0.5*B2_L;
  Real pbr = 0.5*B2_R;
    
  Real hl  = (etot_L + pres_L + pbl)/dens_L;
  Real hr  = (etot_R + pres_R + pbr)/dens_R;
  Real hroe  = (dsl*hl + dsr*hr)*ids;

  Real di  = 1.0/droe;
  Real vsq = sq(uxroe) + sq(uyroe) + sq(uzroe);
  Real btsq = sq(byroe) + sq(bzroe);
  Real bt_startsq = ((kGamma - 1) - (kGamma - 2)*y)*btsq;
  Real vaxsq = Bx*Bx*di;
  Real hp = hroe - (vaxsq + btsq*di);
  Real twid_asq = ((kGamma - 1)*(hp - 0.5*vsq) - (kGamma - 2)*x);
  
  Real ct2  = bt_startsq*di;
  Real tsum = vaxsq + ct2 + twid_asq;
  Real tdif = vaxsq + ct2 - twid_asq;
  Real cf2_cs2 = sqrtR(sq(tdif) + 4.0*twid_asq*ct2);
  Real cfsq = 0.5*(tsum + cf2_cs2);
  Real cf = sqrtR(cfsq);

  Real S_L = uxroe - cf;
  Real S_R = uxroe + cf;

  Real asq = sqrtR(kGamma * pres_L/dens_L);
  vaxsq = sq(Bx)/dens_L;
  ct2   = (sq(By_L) + sq(Bz_L))/dens_L;
  Real qsq = vaxsq + ct2 + asq;
  Real tmp = vaxsq + ct2 - asq;
  cfsq = 0.5*(qsq + sqrtR(sq(tmp) + 4.0*asq*ct2));
  Real cfl = fsig*sqrtR(cfsq);

  asq = sqrtR(kGamma * pres_R/dens_R);
  vaxsq = sq(Bx)/dens_R;
  ct2   = (sq(By_R) + sq(Bz_R))/dens_R;
  qsq = vaxsq + ct2 + asq;
  tmp = vaxsq + ct2 - asq;
  cfsq = 0.5*(qsq + sqrtR(sq(tmp) + 4.0*asq*ct2));
  Real cfr = fsig*sqrtR(cfsq);

  if (velx_R + cfr > S_R) S_R = velx_R + cfr;
  if (velx_L - cfl < S_L) S_L = velx_L - cfl;
#endif
  
  Real pT_L = pres_L + 0.5 * B2_L;
  Real pT_R = pres_R + 0.5 * B2_R;
  
  Real S_M  = (S_R - velx_R)*momx_R - (S_L - velx_L)*momx_L - pT_R + pT_L;
  S_M *= 1.0/((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);

  Real pT_s = (S_R - velx_R)*dens_R*pT_L - (S_L - velx_L)*dens_L*pT_R;
  pT_s += dens_L*dens_R*(S_R - velx_R)*(S_L - velx_L)*(velx_R - velx_L);
  pT_s *= 1.0/((S_R - velx_R)*dens_R - (S_L - velx_L)*dens_L);
  
  Real velx_L_s  = S_M;
  Real velx_L_ss = S_M;
  
  Real velx_R_s  = S_M;
  Real velx_R_ss = S_M;

  Real B2x = Bx*Bx;

  Real dens_L_s = dens_L * (S_L - velx_L)/(S_L - S_M);
  Real vely_L_s = vely_L - Bx*By_L*(S_M - velx_L)/(dens_L*(S_L - velx_L)*(S_L - S_M) - B2x + TINY);
  Real velz_L_s = velz_L - Bx*Bz_L*(S_M - velx_L)/(dens_L*(S_L - velx_L)*(S_L - S_M) - B2x + TINY);
  Real   By_L_s = By_L * (dens_L*sq(S_L - velx_L) - B2x)/(dens_L*(S_L - velx_L)*(S_L - S_M) - B2x + TINY);
  Real   Bz_L_s = Bz_L * (dens_L*sq(S_L - velx_L) - B2x)/(dens_L*(S_L - velx_L)*(S_L - S_M) - B2x + TINY);

  Real dens_R_s = dens_R * (S_R - velx_R)/(S_R - S_M);
  Real vely_R_s = vely_R - Bx*By_R*(S_M - velx_R)/(dens_R*(S_R - velx_R)*(S_R - S_M) - B2x + TINY);
  Real velz_R_s = velz_R - Bx*Bz_R*(S_M - velx_R)/(dens_R*(S_R - velx_R)*(S_R - S_M) - B2x + TINY);
  Real   By_R_s = By_R * (dens_R*sq(S_R - velx_R) - B2x)/(dens_R*(S_R - velx_R)*(S_R - S_M) - B2x + TINY);
  Real   Bz_R_s = Bz_R * (dens_R*sq(S_R - velx_R) - B2x)/(dens_R*(S_R - velx_R)*(S_R - S_M) - B2x + TINY);

  Real   vB_L   = velx_L  *Bx + vely_L  *By_L   + velz_L  *Bz_L;
  Real   vB_L_s = velx_L_s*Bx + vely_L_s*By_L_s + velz_L_s*Bz_L_s;
  Real etot_L_s = ((S_L - velx_L)*etot_L - pT_L*velx_L + pT_s*S_M + Bx*(vB_L - vB_L_s))/(S_L - S_M);

  Real   vB_R   = velx_R  *Bx + vely_R  *By_R   + velz_R  *Bz_R;
  Real   vB_R_s = velx_R_s*Bx + vely_R_s*By_R_s + velz_R_s*Bz_R_s;
  Real etot_R_s = ((S_R - velx_R)*etot_R - pT_R*velx_R + pT_s*S_M + Bx*(vB_R - vB_R_s))/(S_R - S_M);

  Real dens_L_ss = dens_L_s;
  Real    S_L_s  = S_M - absR(Bx)/sqrtR(dens_L_s);

  Real dens_R_ss = dens_R_s;
  Real    S_R_s  = S_M + absR(Bx)/sqrtR(dens_R_s);

  Real f = 1.0/(sqrtR(dens_L_s) + sqrtR(dens_R_s));
  Real vely_ss = f*(sqrtR(dens_L_s)*vely_L_s + sqrtR(dens_R_s)*vely_R_s + (By_R_s - By_L_s)*signBx);
  Real velz_ss = f*(sqrtR(dens_L_s)*velz_L_s + sqrtR(dens_R_s)*velz_R_s + (Bz_R_s - Bz_L_s)*signBx);

  Real By_ss = f*(sqrtR(dens_L_s)*By_R_s + sqrtR(dens_R_s)*By_L_s + sqrtR(dens_L_s*dens_R_s)*(vely_R_s - vely_L_s)*signBx);
  Real Bz_ss = f*(sqrtR(dens_L_s)*Bz_R_s + sqrtR(dens_R_s)*Bz_L_s + sqrtR(dens_L_s*dens_R_s)*(velz_R_s - velz_L_s)*signBx);
  
  Real vely_L_ss = vely_ss;
  Real velz_L_ss = velz_ss;
  Real   By_L_ss = By_ss;
  Real   Bz_L_ss = Bz_ss;

  Real vely_R_ss = vely_ss;
  Real velz_R_ss = velz_ss;
  Real   By_R_ss = By_ss;
  Real   Bz_R_ss = Bz_ss;

  Real vB_L_ss   = velx_L_ss*Bx + vely_L_ss*By_L_ss + velz_L_ss*Bz_L_ss;
  Real etot_L_ss = etot_L_s - sqrtR(dens_L_s)*(vB_L_s - vB_L_ss)*signBx;

  Real vB_R_ss   = velx_R_ss*Bx + vely_R_ss*By_R_ss + velz_R_ss*Bz_R_ss;
  Real etot_R_ss = etot_R_s + sqrtR(dens_R_s)*(vB_R_s - vB_R_ss)*signBx;

  Real Fdens_L = dens_L*velx_L;
  Real Fmomx_L = momx_L*velx_L + pT_L        - B2x;
  Real Fmomy_L = momy_L*velx_L               - Bx*By_L;
  Real Fmomz_L = momz_L*velx_L               - Bx*Bz_L;
  Real Fetot_L = etot_L*velx_L + pT_L*velx_L - Bx*vB_L; 
  Real Fby_L   = By_L  *velx_L               - Bx*vely_L;
  Real Fbz_L   = Bz_L  *velx_L               - Bx*velz_L;

  Real Fdens_R = dens_R*velx_R;
  Real Fmomx_R = momx_R*velx_R + pT_R        - B2x;
  Real Fmomy_R = momy_R*velx_R               - Bx*By_R;
  Real Fmomz_R = momz_R*velx_R               - Bx*Bz_R;
  Real Fetot_R = etot_R*velx_R + pT_R*velx_R - Bx*vB_R;
  Real Fby_R   = By_R  *velx_R               - Bx*vely_R;
  Real Fbz_R   = Bz_R  *velx_R               - Bx*velz_R;

  Real momx_L_s  = dens_L_s *velx_L_s;
  Real momy_L_s  = dens_L_s *vely_L_s;
  Real momz_L_s  = dens_L_s *velz_L_s;
  
  Real momx_L_ss = dens_L_ss*velx_L_ss;
  Real momy_L_ss = dens_L_ss*vely_L_ss;
  Real momz_L_ss = dens_L_ss*velz_L_ss;

  Real momx_R_s  = dens_R_s *velx_R_s;
  Real momy_R_s  = dens_R_s *vely_R_s;
  Real momz_R_s  = dens_R_s *velz_R_s;
  
  Real momx_R_ss = dens_R_ss*velx_R_ss;
  Real momy_R_ss = dens_R_ss*vely_R_ss;
  Real momz_R_ss = dens_R_ss*velz_R_ss;


  if (S_L > 0) {
    Fdens = Fdens_L;
    Fmomx = Fmomx_L;
    Fmomy = Fmomy_L;
    Fmomz = Fmomz_L;
    Fetot = Fetot_L;
    F_By  = Fby_L;
    F_Bz  = Fbz_L;

  } else if (S_L <= 0 && 0 <= S_L_s) {
    Fdens = Fdens_L + S_L*dens_L_s - S_L*dens_L;
    Fmomx = Fmomx_L + S_L*momx_L_s - S_L*momx_L;
    Fmomy = Fmomy_L + S_L*momy_L_s - S_L*momy_L;
    Fmomz = Fmomz_L + S_L*momz_L_s - S_L*momz_L;
    Fetot = Fetot_L + S_L*etot_L_s - S_L*etot_L;
    F_By  = Fby_L   + S_L*By_L_s   - S_L*By_L;
    F_Bz  = Fbz_L   + S_L*Bz_L_s   - S_L*Bz_L;

  } else if (S_L_s <= 0 && 0 <= S_M) {
    Fdens = Fdens_L + S_L_s*dens_L_ss - (S_L_s - S_L)*dens_L_s - S_L*dens_L;
    Fmomx = Fmomx_L + S_L_s*momx_L_ss - (S_L_s - S_L)*momx_L_s - S_L*momx_L;
    Fmomy = Fmomy_L + S_L_s*momy_L_ss - (S_L_s - S_L)*momy_L_s - S_L*momy_L;
    Fmomz = Fmomz_L + S_L_s*momz_L_ss - (S_L_s - S_L)*momz_L_s - S_L*momz_L;
    Fetot = Fetot_L + S_L_s*etot_L_ss - (S_L_s - S_L)*etot_L_s - S_L*etot_L;
    F_By  = Fby_L   + S_L_s*By_L_ss   - (S_L_s - S_L)*By_L_s   - S_L*By_L;
    F_Bz  = Fbz_L   + S_L_s*Bz_L_ss   - (S_L_s - S_L)*Bz_L_s   - S_L*Bz_L;

  } else if (S_M <= 0 && 0 <= S_R_s) {
    Fdens = Fdens_R + S_R_s*dens_R_ss - (S_R_s - S_R)*dens_R_s - S_R*dens_R;
    Fmomx = Fmomx_R + S_R_s*momx_R_ss - (S_R_s - S_R)*momx_R_s - S_R*momx_R;
    Fmomy = Fmomy_R + S_R_s*momy_R_ss - (S_R_s - S_R)*momy_R_s - S_R*momy_R;
    Fmomz = Fmomz_R + S_R_s*momz_R_ss - (S_R_s - S_R)*momz_R_s - S_R*momz_R;
    Fetot = Fetot_R + S_R_s*etot_R_ss - (S_R_s - S_R)*etot_R_s - S_R*etot_R;
    F_By  = Fby_R   + S_R_s*By_R_ss   - (S_R_s - S_R)*By_R_s   - S_R*By_R;
    F_Bz  = Fbz_R   + S_R_s*Bz_R_ss   - (S_R_s - S_R)*Bz_R_s   - S_R*Bz_R;

  } else if (S_R_s <= 0 && 0 <= S_R) {
    Fdens = Fdens_R + S_R*dens_R_s - S_R*dens_R;
    Fmomx = Fmomx_R + S_R*momx_R_s - S_R*momx_R;
    Fmomy = Fmomy_R + S_R*momy_R_s - S_R*momy_R;
    Fmomz = Fmomz_R + S_R*momz_R_s - S_R*momz_R;
    Fetot = Fetot_R + S_R*etot_R_s - S_R*etot_R;
    F_By  = Fby_R   + S_R*By_R_s   - S_R*By_R;
    F_Bz  = Fbz_R   + S_R*Bz_R_s   - S_R*Bz_R;

  } else {
    Fdens = Fdens_R;
    Fmomx = Fmomx_R;
    Fmomy = Fmomy_R;
    Fmomz = Fmomz_R;
    Fetot = Fetot_R;
    F_By  = Fby_R;
    F_Bz  = Fbz_R;

  }

}

