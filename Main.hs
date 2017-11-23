-- file: Astro/Main.hs
-- For running benchmarks

module Main (main)
	where

import Astro.Benchmarks

main = do
	keplerPropagateBM 1000
