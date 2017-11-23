-- file: Astro/Leapfrog.hs
-- Collection of abstract leapfrog integrator utilities.

module Astro.Leapfrog where

--type StepF = Double->a->a

lf2_qpq :: Double->(Double->a->a)->(Double->a->a)->a->a
lf2_qpq h qstep pstep = qstep (h/2) . pstep h . qstep (h/2)

lf2_pqp :: Double->(Double->a->a)->(Double->a->a)->a->a
lf2_pqp h qstep pstep = pstep (h/2) . pstep h . qstep (h/2)
