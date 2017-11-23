{-# LANGUAGE BangPatterns #-}
-- file: Astro/Misc.hs
-- Miscellaneous utilities

module Astro.Misc where

import Data.List
import Debug.Trace
import Numeric.LinearAlgebra
import Numeric (showEFloat)
import Text.Printf

{-
 - General math
 -}

type Polynomial = [Double]

-- Evaluate real polynomial given as list of coefficients [a0,a1,...]
evalPoly :: Polynomial->Double->Double
evalPoly p x = sum $ zipWith (*) p $ map (x^) [0..]

-- Derivative of a (real) polynomial
diffPoly :: Polynomial->Polynomial
diffPoly [] = []
diffPoly (a0:[]) = []
diffPoly (a0:as) = zipWith (*) as [1..]

-- Degree of a polynomial
degPoly :: Polynomial->Int
degPoly [] = 0
degPoly (a:[]) = if a /= 0 then 1 else 0
degPoly p = if l == [] then 0 else last l
	where
	l = Data.List.findIndices (/= 0) p

-- Laguerre polynomial root finder
-- Returns Nothing if fails to converge to a real root within maxIt
rootLaguerre (tol,maxIt) p x0  = result
	where
	result = helper x0 0
	p' = diffPoly p
	p'' = diffPoly p'
	deg = fromIntegral $ degPoly p
	helper x n
		| isNaN x = Nothing
		| abs (evalPoly p x) <= tol = Just x
		| n >= maxIt = Nothing -- `debug` ("maxit x,nx:"++show(x,nx))
		| otherwise = helper nx (n+1) 
                {-
                        `debug` (printf ("g:n %g h: %g sqa: %g sq: %g x:"
                                ++"%g y: %g") g h sqarg sq nx y)
                                -}
		where
		g = evalPoly p' x / evalPoly p x
		h = g^2 - evalPoly p'' x / evalPoly p x
		sqarg = (deg-1)*(deg*h-g^2)
		sq = sqrt sqarg
		denom = if g>0 then g+sq else g-sq
		nx = x - deg/denom -- `debug` ("n,denom:"++show(n,denom))
                y = evalPoly p nx

{-
 - Vector Math
 -}

-- Standard vector norms
norm :: Vector Double -> Double
--norm v = pnorm PNorm2 v
-- XXX: For 3-vectors this is quite a bit faster
norm !v = sqrt $ v@>0*v@>0 + v@>1*v@>1 + v@>2*v@>2 

norm3D :: Vector Double -> Double
norm3D = sqrt . sqrNorm3D

sqrNorm3D :: Vector Double -> Double
sqrNorm3D v = v@>0*v@>0 + v@>1*v@>1 + v@>2*v@>2 

-- Normalize vector
unitV v = scale (recip $ norm v) v

-- 3-vector cross product
cross :: Vector Double->Vector Double->Vector Double
cross !a !b =  3|>[a@>1 * b@>2 - a@>2 * b@>1 
		,-a@>0 * b@>2 + a@>2 * b@>0 
		,a@>0 * b@>1 - a@>1 * b@>0] 

-- Rodrigues' rotation formula
rodr :: Vector Double->Double->Vector Double->Vector Double
rodr z th vec = scale (cos th) vec + scale (sin th) (cross z vec)
		+ scale ((z<.>vec)*(1-cos th)) z

{-
 - List and tuple manipulations
 -}

thrd (_,_,a) = a
frth (_,_,_,a) = a

-- Give all pairings, excluding (p,p) and when (p1,p2) excluding (p2,p1)
-- Courtesy of Jussi Knuuttila
chooseTwo :: [a] -> [(a,a)]
chooseTwo as = concatMap pairs (tails as)
    where pairs [] = []
          pairs (a:as) = zip (repeat a) as

-- Take every nth of list
takenth 0 _ = []
takenth _ [] = []
takenth n list = 
	iter list []
	where
	iter [] acc = acc
	iter list acc =
		iter (tail newList) (acc ++ [head newList])
		where newList = drop (n-1) list

-- Debugging and printing functions
debug = flip trace
showFloat f = showEFloat (Just 4) f ""
