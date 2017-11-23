-- file: Astro/RunTests.hs
-- Run QuickCheck tests

import Astro.Stumpff
import Astro.Tests
import Test.QuickCheck

main = do
	putStrLn "stumpff"
	quickCheckResult prop_stumpff_c_0
	quickCheckResult prop_stumpff_c_rec
	quickCheckResult prop_stumpff_G_0
	quickCheckResult prop_stumpff_G_rec
	putStrLn "orbital elements"
	quickCheckResult prop_oe_rv_id
		


