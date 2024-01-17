/*
 *  spectrum.cpp
 *  
 *
 *  Created by tmo22 on 2/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectrum.h"

Doub square(Int j,Int n) {return 1.;}

Doub bartlett(Int j,Int n) {return 1.-abs(2.*j/(n-1.)-1.);}

Doub welch(Int j,Int n) {return 1.-SQR(2.*j/(n-1.)-1.);}

void Slepian::filltable() {
	const Doub EPS = 1.e-10, PI = 4.*atan(1.);
	Doub xx,xnew,xold,sw,ppp,ddd,sum,bet,ssub,ssup,*u;
	Int i,j,k,nl;
	VecDoub dg(m2),dgg(m2),gam(m2),sup(m2-1),sub(m2-1);
	sw = 2.*SQR(sin(jres*PI/m2));
	dg[0] = 0.25*(2*m2+sw*SQR(m2-1.)-1.);
	for (i=1;i<m2;i++) {
		dg[i] = 0.25*(sw*SQR(m2-1.-2*i)+(2*(m2-i)-1.)*(2*i+1.));
		sub[i-1] = sup[i-1] = -i*(Doub)(m2-i)/2.;
	}
	xx = -0.10859 - 0.068762/jres + 1.5692*jres;
	xold = xx + 0.47276 + 0.20273/jres - 3.1387*jres;
	for (k=0; k<kt; k++) {
		u = &dpss[k][0];
		for (i=0;i<20;i++) {
			pp = 1.;
			p = dg[0] - xx;
			dd = 0.;
			d = -1.;
			for (j=1; j<m2; j++) {
				ppp = pp; pp = p;
				ddd = dd; dd = d;
				p = pp*(dg[j]-xx) - ppp*SQR(sup[j-1]);
				d = -pp + dd*(dg[j]-xx) - ddd*SQR(sup[j-1]);
				if (abs(p)>1.e30) renorm(-100);
				else if (abs(p)<1.e-30) renorm(100);
			}
			xnew = xx - p/d;
			if (abs(xx-xnew) < EPS*abs(xnew)) break;
			xx = xnew;
		}
		xx = xnew - (xold - xnew);
		xold = xnew;
		for (i=0;i<m2;i++) dgg[i] = dg[i] - xnew;
		nl = m2/3;
		dgg[nl] = 1.;
		ssup = sup[nl]; ssub = sub[nl-1];
		u[0] = sup[nl] = sub[nl-1] = 0.;
		bet = dgg[0];
		for (i=1; i<m2; i++) {
			gam[i] = sup[i-1]/bet;
			bet = dgg[i] - sub[i-1]*gam[i];
			u[i] = ((i==nl? 1. : 0.) - sub[i-1]*u[i-1])/bet;
		}
		for (i=m2-2; i>=0; i--) u[i] -= gam[i+1]*u[i+1];
		sup[nl] = ssup; sub[nl-1] = ssub;
		sum = 0.;
		for (i=0; i<m2; i++) sum += SQR(u[i]);
		sum = (u[3] > 0.)? sqrt(sum) : -sqrt(sum);
		for (i=0; i<m2; i++) u[i] /= sum;
	}
}
