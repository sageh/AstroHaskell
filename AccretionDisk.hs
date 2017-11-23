-- file:Astro/AccretionDisk.hs
-- Numerical tools for accretion disk analysis

module Astro.AccretionDisk
	where

import Control.Monad
import System.IO.Unsafe
import System.Random
import Numeric.LinearAlgebra

import Astro.Misc

seed = 1337

rndDoubles :: Int -> (Double,Double) -> IO [Double]
rndDoubles n (min,max) = replicateM n $ randomRIO (min,max)

rndDouble :: (Double, Double) -> IO Double
rndDouble (min, max) = randomRIO (min,max) :: IO Double

-- Random position and corresponding velocity in a flat xy-plane Keplerian disk
-- with CCW rotation when viewed towards negative z-axis
rndDiskPosVel :: (Double, Double) -> Double -> IO (Vector Double, Vector Double)
rndDiskPosVel (minR,maxR) gpar = 
	do
	r <- rndDouble (minR,maxR)
	th <- rndDouble (0, 2*pi)
	let rv = fromList [ r * cos th, r * sin th, 0]
	    vv = fromList [-rv@>1, rv@>0, 0.0]
	    vv' = scale (sqrt gpar/norm rv) vv
	return (rv, vv')

-- Create a flat disk (planar), rotating in a CCW sense when viewed
-- towards negative z-axis
flatDisk :: (Double, Double) -> Double -> Int 
		-> [(Vector Double, Vector Double)]
flatDisk (minR, maxR) gpar nparts = unsafePerformIO $
	do 
		setStdGen $ mkStdGen seed
		replicateM nparts $ rndDiskPosVel (minR,maxR) gpar

	
