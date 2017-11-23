-- file: Astro/Benchmarks.hs
-- Benchmarks for some of the functions

module Astro.Benchmarks
	where

import Astro.Stumpff
import Astro.Kepler
import Astro.Elements

import Numeric.GSL
import Numeric.LinearAlgebra
import Control.Monad

keplerPropagateBM :: Int -> IO ()
keplerPropagateBM nsteps = step 0 (rv,vv) 
	where
	(rv, vv) 	= (3|>[1,0,0], 3|>[0,0.1,0])
	oe 		= posvelToOrbel rv vv
	h  		= 10*2*pi/(fromIntegral nsteps::Double)
	step n (r,v) = do
		putStrLn $ (show r)
		(when (n < nsteps) $ step (n+1) newS)
		where 
		newS = keplerPropagate r v h (1e-15,100)
