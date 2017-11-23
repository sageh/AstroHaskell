-- file: Astro/Elements.hs
-- Orbital element calculations.

module Astro.Elements
	where

import Astro.Misc
import Astro.Kepler
import Numeric.LinearAlgebra

-- Standard orbital elements from position and velocity, with given
-- gravitational parameter
posvelToOrbelGp :: Vector Double -> Vector Double -> Double -> [Double]
	--(Monad m) => Vector Double -> Vector Double -> Double -> m [Double]
posvelToOrbelGp rv vv' gpar = 
	{-
	if (vv'<.>vv')/2 - gpar / (norm rv) >= 0 then
		fail "Not an elliptic orbit"
	else
	-}
	let 	vv 	= scalar (1/sqrt gpar) * vv'
		r 	= norm rv
		alpha 	= 2/r - vv<.>vv
		eta 	= rv <.> vv
		zeta	= 1 - alpha * r
		hv	= cross rv vv
		h	= norm hv 
		ev	= cross vv hv - scale (1/r) rv
		b	= sqrt $ (hv@>0)^2 + (hv@>1)^2
		m0	= atan2 (eta * sqrt alpha) zeta - eta * sqrt alpha
		inc	= atan2 b (hv@>2)
		ascnode = atan2 (hv@>0) (-hv@>1)
		argperi = atan2 ((ev@>2)*h) (ev@>1*hv@>0-ev@>0*hv@>1)
	in	[1/alpha, norm ev, inc, ascnode, argperi, m0]
		--return $ [1/alpha, norm ev, inc, ascnode, argperi, m0]
posvelToOrbel rv vv = posvelToOrbelGp rv vv 1

-- Position and velocity from orbital elements, with given gravitational
-- parameter
orbelToPosvelGp :: [Double] -> Double -> (Vector Double, Vector Double)
	--(Monad m) => [Double] -> Double -> m (Vector Double, Vector Double)
orbelToPosvelGp oe gpar =
	{-
	if (length oe < 6) then
		fail "Need 6 orbital elements"
	else
	-}
	let 	[a, ecc, inc, o, w, m0] = oe
		pv	= fromList [ cos w *cos o - sin w *sin o *cos inc 
		 	    	   , cos w *sin o + sin w *cos o *cos inc 
		 	    	   , sin w *sin inc ]
		qv	= fromList [ -sin w *cos o  - cos w *sin o *cos inc 
		 	   	   , -sin w *sin o  + cos w *cos o *cos inc 
		 	    	   ,  cos w *sin inc ]
		wv	= fromList [  sin o *sin inc 
				   , -cos o *sin inc 
				   ,  cos inc  ]
		ea	= solveEKepler ecc m0
		av	= scalar a * pv
		bv	= scalar (sqrt (1-ecc^2)) * cross wv av
		rv	= scalar (cos ea - ecc) * av 
					+ scalar (sin ea) * bv 
		vv	= scalar ((sqrt gpar)/(norm rv * sqrt a)) * 
				( scalar (sin ea) * (-av) 
					+ scalar (cos ea) * bv )
	in
		(rv, vv)
		--return $ (rv, vv)
orbelToPosvel oe = orbelToPosvelGp oe 1
