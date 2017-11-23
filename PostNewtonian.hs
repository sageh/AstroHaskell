-- file: Astro/PostNewtonian.hs
-- Post-Newtonian approximations for two-body problem

module Astro.PostNewtonian
	where

import Astro.Misc
import Numeric.LinearAlgebra

import Data.List (foldl')

-- Calculate Post-Newtonian acceleration functions for non-spin
-- dependent terms
mkNonSpinPNs :: Double->Double->Double->Double->
	[(Vector Double->Vector Double->Vector Double)]
mkNonSpinPNs gc c m1 m2 = 	[ a1
				, a2
				, a2_5
				] 
			where
	-- Keplerian contribution
	--a0 :: Vector Double -> Vector Double -> Vector Double
	--a0 rv vv = scale (-gm / r^3) rv where r = norm rv

	a1 :: Vector Double -> Vector Double -> Vector Double
	a1 rv vv = {-# SCC "a1" #-} 
		scale (-gmperc2 / r^2) $
			scale nbrac nv + scale vbrac vv
		where
			nbrac = -2*(2+eta)*gm/r + (1+3*eta)*v^2 - 3/2*eta*rd^2 
			vbrac = -2*(2-eta)*rd
			r	= norm rv
			v	= norm vv
			nv	= unitV rv
			rd	= (rv <.> vv)/r

	a2 :: Vector Double -> Vector Double -> Vector Double
	a2 rv vv = {-# SCC "a2" #-} 
		scale (- gmperc2 / (c^2 * r^2)) $
			scale (nbrac/r) rv + scale ((-rd/2) * vbrac) vv
		where
			nbrac = 3/4*(12+29*eta)*(gm/r)^2 + eta*(3-4*eta)*v^4
				+ 15/8*eta*(1-3*eta)*rd^4
				- 3/2*eta*(3-4*eta)*v^2*rd^2
				- 1/2*eta*(13-4*eta)*(gm/r)*v^2
				- (2+25*eta+2*eta^2)*(gm/r)*rd^2
			vbrac = eta*(15+4*eta)*v^2 - (eta+41*eta+8*eta^2)*(gm/r)
				- 3*eta*(3+2*eta)*rd^2
			r	= norm rv
			v	= norm vv
			rd	= (rv <.> vv)/r

	a2_5 :: Vector Double -> Vector Double -> Vector Double
	a2_5 rv vv = {-# SCC "a2_5" #-} 
		scale (8/15*(gmperc2)^2*eta/(c*r^2)) $
			scale (nbrac*rd/r) rv + scale (-vbrac) vv
		where
			nbrac 	= 9*v^2 + 17*gm/r
			vbrac 	= 3*v^2 + 9*gm/r
			r	= norm rv
			v	= norm vv
			rd	= (rv <.> vv)/r
	gm 	= gc * m
	gmperc2 = gm / c ^2
	m 	= m1+m2
	eta 	= m1*m2/m^2

-- Calculate Post-Newtonian acceleration functions for spin
-- dependent terms
--mkSpinPNs :: Double->Double->Double->Double->Double->Double->
--	[(Vector Double->Vector Double->Vector Double->Vector Double)]
mkSpinPNs gc c m1 m2 chi q = 	[ spinOrbit
				, quadMonoInteraction
				] where
	--spinOrbit :: Vector Double -> Vector Double -> Vector Double
	--		-> Vector Double
	spinOrbit rv vv s = {-# SCC "SO" #-} 
		scale (gm^2/(r^3*c^3)*(1+sqrt (1-4*eta))/4*chi) $
				scale nbrac nv + scale (nsbrac) (cross nv s')
					+ scale (-vsbrac) (cross vv s')
		where
			nbrac = 12 * s' <.> (cross nv vv)
			nsbrac = (9 + 3 * sqrt (1-4*eta))*rd
			vsbrac = 7 + sqrt(1-4*eta)
			r	= norm rv
			v	= norm vv
			nv	= unitV rv
			rd	= (rv <.> vv)/r
			s'	= spinCylToSpher s

	--quadMonoInteraction :: Vector Double -> Vector Double -> Vector Double
	--			-> Vector Double
	quadMonoInteraction rv vv s = {-# SCC "Q" #-}
		-- FIXME: Minus or no minus? Depends on whether
		-- coordinates fixed on m1 or m2, see OJ287 paper, minus
		-- in paper.
		scale (-q*chi^2*3*gc^3*m1^2*m/(2*c^4*r^4)) $
			scale nbrac nv + scale (-2 * nv <.> s') s'
		where
			nbrac = 5 * (nv <.> s')^2 - 1
			r	= norm rv
			v	= norm vv
			nv	= unitV rv
			s'	= spinCylToSpher s
	gm 	= gc * m
	m 	= m1+m2
	eta 	= m1*m2/m^2

-- Time derivative of spin unit vector, in spherical coordinates
spinDt :: Double->Double->Double->Double->Double->Double->
	Vector Double->Vector Double->Vector Double->Vector Double
spinDt gc c m1 m2 chi q rv vv s = cross omega s 
	where
	omega = scale (gm*eta/(2*c^2*r^2)*(7+sqrt (1-4*eta))
				/(1+sqrt (1-4*eta))) (cross nv vv)
	r = norm rv
	nv = scale (1/r) rv
	gm 	= gc * m
	m 	= m1+m2
	eta 	= m1*m2/m^2

-- Time derivatives of spin unit vector components in cylindrical
-- coordinates, with x as the axisymmery axis
spinDtCyl :: Double->Double->Double->Double->Double->Double->
	Vector Double->Vector Double->[(Double,Double)->Double]
spinDtCyl gc c m1 m2 chi q rv vv = [dtheta, dxi]
	where
	dtheta (th,xi) = 
		-xiPerRho*(dhs@>1*cos th + dhs@>2*sin th) + dhs@>0
		where
		xiPerRho = if abs xi >= 1 then 0 else xi/sqrt (1-xi^2)
	dxi (th,xi) =
		rho * (dhs@>1*sin th - dhs@>2 * cos th)
		where
		rho = if abs xi >= 1 then 0 else sqrt $ 1 - xi^2
	r = norm rv
	nv = scale (1/r) rv
	dhs = omega
	omega = scale (gm*eta/(2*c^2*r^2)*(7+sqrt (1-4*eta))
				/(1+sqrt (1-4*eta))) (cross nv vv)
	gm 	= gc * m
	m 	= m1+m2
	eta 	= m1*m2/m^2


-- Convert spin from cylindrical to cartesian representation
-- Axisymmetry axis is x-axis
spinCylToSpher :: (Double,Double)->Vector Double
spinCylToSpher (theta,xi) = --3|>[rho * cos theta, rho*sin theta, xi]
	3|>[xi, rho * cos theta, rho*sin theta]
	where rho = if abs xi >= 1 then 0 else sqrt $ 1 - xi^2

-- Post-Newtonian acceleration functions for highest supported degree
nonSpinPNAcc gc c m1 m2 rv vv = 
	sum [f rv vv | f <- mkNonSpinPNs gc c m1 m2]
spinPNAcc gc c m1 m2 chi q rv vv s = 
	sum [f rv vv s | f <- mkSpinPNs gc c m1 m2 chi q]

totalPNAcc gc c m1 m2 chi q rv vv s = ((f1 rv vv) + (f2 rv vv s))
	where (f1, f2) = (nonSpinPNAcc gc c m1 m2, 
				spinPNAcc gc c m1 m2 chi q)

-- Explicitly calculated total acceleration
totalPNAccExp gc c m1 m2 chi q rv vv s = 
	a1 + spinOrbit + quadMono + a2 + a2_5
	where
	a1 = {-# SCC "a1" #-} 
		scale (-gmperc2 / r^2) $
			scale nbrac nv + scale vbrac vv
		where
			nbrac = -2*(2+eta)*gm/r + (1+3*eta)*v^2 - 3/2*eta*rd^2 
			vbrac = -2*(2-eta)*rd

	a2 = {-# SCC "a2" #-} 
		scale (- gmperc2 / (c^2 * r^2)) $
			scale nbrac nv + scale ((-rd/2) * vbrac) vv
		where
			nbrac = 3/4*(12+29*eta)*(gm/r)^2 + eta*(3-4*eta)*v^4
				+ 15/8*eta*(1-3*eta)*rd^4
				- 3/2*eta*(3-4*eta)*v^2*rd^2
				- 1/2*eta*(13-4*eta)*(gm/r)*v^2
				- (2+25*eta+2*eta^2)*(gm/r)*rd^2
			vbrac = eta*(15+4*eta)*v^2 - (eta+41*eta+8*eta^2)*(gm/r)
				- 3*eta*(3+2*eta)*rd^2

	a2_5 = {-# SCC "a2_5" #-} 
		scale (8/15*(gmperc2)^2*eta/(c*r^2)) $
			scale (nbrac*rd) nv + scale (-vbrac) vv
		where
			nbrac 	= 9*v^2 + 17*gm/r
			vbrac 	= 3*v^2 + 9*gm/r

	spinOrbit = {-# SCC "SO" #-} 
		scale (gm^2/(r^3*c^3)*(1+sqrt (1-4*eta))/4*chi) $
				scale nbrac nv + scale (nsbrac) (cross nv s')
					+ scale (-vsbrac) (cross vv s')
		where
			nbrac = 12 * s' <.> (cross nv vv)
			nsbrac = (9 + 3 * sqrt (1-4*eta))*rd
			vsbrac = 7 + sqrt(1-4*eta)

	quadMono = {-# SCC "Q" #-}
		scale (-q*chi^2*3*gc^3*m1^2*m/(2*c^4*r^4)) $
			scale nbrac nv + scale (-2 * nv <.> s') s'
		where
			nbrac = 5 * (nv <.> s')^2 - 1
	r	= norm rv
	v	= norm vv
	s'	= spinCylToSpher s
	nv	= scale (1/r) rv
	rd	= (rv <.> vv)/r
	gm 	= gc * m
	gmperc2 = gm / c ^2
	m 	= m1+m2
	eta 	= m1*m2/m^2
