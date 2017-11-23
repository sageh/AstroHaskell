-- file: Astro/Stumpff.hs
-- Stumpff c & G functions

module Astro.Stumpff ( stumpff_c
		     , stumpff_G
		     )
	where

import Complex

-- Helper functions
factlist = 1 : (zipWith (*) [1..] factlist)
fac :: Int -> Int
fac n 
	| n >= 0 = factlist !! n
	| n < 0  = 0 -- to facilitate negative k (for whatever reason)

-- Calculate Stumpff c function of degree k, recursing for k > 4
-- FIXME: Numerical instabilities when k -> inf
stumpff_c_cplx :: Int -> Double -> Double
stumpff_c_cplx _ 0.0 =  0.0
stumpff_c_cplx k x 
	| k == 0 = realPart $ cos (sqrt x')
	| k == 1 = realPart $  sin (sqrt x') / sqrt x'
	| k == 2 = realPart $ (1.0 - cos (sqrt x')) / x'
	| k == 3 = realPart $  (sqrt x' - sin (sqrt x')) / (x' * sqrt x')
	| otherwise = (1/fromIntegral (fac (k-2)) - stumpff_c_cplx (k-2) x)/x
	where
		x' :: Complex Double
		x' =  realToFrac x

-- Calculate Stumpff G function of degree k with parameter beta,
-- recursing for k > 4
stumpff_G_cplx :: Int -> Double -> Double -> Double
stumpff_G_cplx _ _ 0.0 = 0
stumpff_G_cplx k beta x = 
	(x ** fromIntegral k) * stumpff_c_cplx k (beta * (x^2))

stumpff_c = stumpff_c_real
stumpff_G = stumpff_G_real

-- Same without complex numbers
stumpff_c_real :: Int -> Double -> Double
stumpff_c_real _ 0.0 = 0.0
stumpff_c_real 0 x = if x>=0 
	then cos $ sqrt x else cos $ sqrt x
stumpff_c_real 1 x = if x>=0
	then sin (sqrt x) / sqrt x else sinh (sqrt x) / sqrt x
stumpff_c_real 2 x = if x>=0 
	then (1 - cos (sqrt x)) / x else (-1 + cosh (sqrt x))/x
stumpff_c_real 3 x = if x>=0 
		then (sqrt x - sin (sqrt x)) / x**(3/2)
		else (-sqrt x + sinh (sqrt x))/ x**(3/2)
stumpff_c_real k x =
	(1/fromIntegral (fac (k-2)) - stumpff_c_real (k-2) x)/x

stumpff_G_real :: Int -> Double -> Double -> Double
stumpff_G_real _ _ 0.0 = 0
stumpff_G_real k beta x = 
	(x ** fromIntegral k) * stumpff_c_real k (beta * (x^2))
