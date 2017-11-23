-- file: Astro/Tests.hs
-- Tests for some of the library functions

module Astro.Tests
	where

import Astro.Stumpff
import Astro.Kepler
import Astro.Elements

import Control.Monad
import Numeric.LinearAlgebra
import Test.QuickCheck
import Test.QuickCheck.Arbitrary

-- Stumpff c functions
stumpff_c_rec_lhs k x = x * stumpff_c k x
stumpff_c_rec_rhs k x = 1.0/fromIntegral (fac (k-2)) - stumpff_c (k-2) x
	where fac n = product [1..n]
prop_stumpff_c_0 k = stumpff_c k 0.0 == 0
prop_stumpff_c_rec k x = 
	let lhs = stumpff_c_rec_lhs k x
	    rhs = stumpff_c_rec_rhs k x
	in
	(k >= 2 && k <= 35 && x /= 0) ==> abs $ (rhs-lhs)/lhs <= 1e-7
	where fac n = product [1..n]
prop_stumpff_G_0 k beta = (stumpff_G k beta 0.0 == 0)
prop_stumpff_G_rec k beta x = 
	(k >= 2) ==> (x * stumpff_G k beta x == 
		x^(k-2)/fromIntegral (fac (k-2)) - stumpff_G (k-2) beta x)
	where fac n = product [1..n]

-- Orbital elements
{-
instance Arbitrary (Vector t) where
	arbitrary = liftM fromList (vectorOf 3 arbitrarySizedFractional)
-}
--prop_oe_rv_id :: Vector Double -> Vector Double -> Double -> Property
prop_oe_rv_id :: [Double] -> [Double] -> Double -> Property
prop_oe_rv_id r v gpar = 
	property $ (length r == 3 && length v == 3) ==> (rv == rv' && vv == vv')
	where
		rv = fromList r
		vv = fromList v
		(rv',vv') = orbelToPosvel_gp (posvelToOrbel_gp rv vv gpar) gpar
	
