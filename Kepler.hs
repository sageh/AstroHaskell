-- file: Astro/Kepler.hs
-- Various forms of Kepler's equation with efficient solvers

module Astro.Kepler ( eKepler
		    , gKepler
		    , solveEKepler
		    , solveGKepler
		    , keplerPropagate
		    )
	where

import Astro.Stumpff
import Astro.Misc

import Complex
import Debug.Trace
import Numeric.GSL
import Numeric.LinearAlgebra
import Text.Printf (printf)

-- Shorthands and such
type CD = Complex Double
type VD = Vector Double

-- Elliptic Kepler's equation
eKepler :: (Floating a) => a -> a -> a
eKepler e ea = ea - e * sin ea

-- General Variable Formulation of Kepler's eq
gKepler :: Double -> Double -> Double -> Double -> Double -> Double
gKepler r0 eta0 zeta0 beta x = 
	r0 * x + eta0 * stumpff_G 2 beta x + zeta0 * stumpff_G 3 beta x

-- Solve Elliptic Kepler's equation iteratively
solveEKepler :: (Floating a) => a -> a -> a
solveEKepler e m =
	iterate 0 m m
	where 
	--iterate n ea m | trace ("n,E: " ++ show n ++ " " ++ show ea) False = undefined
	iterate n ea m 
		| ea' == ea || n > maxit = ea'
		| otherwise = iterate (n+1) ea' m
		where
		    fd = 1.0 - e * cos ea
		    fdd = e * sin ea
		    f = eKepler e ea - m
		    ea' = ea - f / sqrt (fd^2 - f * fdd)
	maxit = 10

-- Solve General Kepler's eq iteratively
--solveGKepler :: Double -> Double -> Double -> Double -> Double -> Double -> Double
solveGKepler r0 eta0 zeta0 beta x0 t (tol,maxIt) = 
	if abs lastErr > tol
		then head sol `debug` 
			("solveGKepler: iter. incomp., error: "
				++(show lastErr))
		else head sol
	where
		(sol,path) = root Hybrids tol maxIt f [x0]
		lastErr = (last.last) $ toLists path
		f vx = [gKepler r0 eta0 zeta0 beta (head vx) - t]

-- Solve r, v given r0, v0, gravitational parameter, delta_t
--keplerPropagate ::  VD -> VD -> Double -> (Double, Int) -> (VD,VD)
keplerPropagate dt (tol,maxIt) rv vv =
	(scale f rv + scale g vv , scale fd rv + scale gd vv)
	where (r0, v0) = (norm rv, norm vv) 
	      tau = dt -- working in units where G(m_1+m_2) = 1
	      eta0 = rv <.> vv
	      beta = 2/r0 - v0^2
	      zeta0 = 1 - beta * r0
	      x = solveGKepler r0 eta0 zeta0 beta 0.0 tau (tol,maxIt)
	      r = r0 + eta0 * stumpff_G 1 beta x + zeta0 * stumpff_G 2 beta x
	      f = 1 - stumpff_G 2 beta x / r0
	      g = tau - stumpff_G 3 beta x
	      fd = - stumpff_G 1 beta x / (r * r0)
	      gd = 1 - stumpff_G 2 beta x / r
