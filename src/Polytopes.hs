module Polytopes where

import Quaternion

class Phi a where
  phi :: a

instance Phi Float where
  phi = (sqrt 5 + 1) / 2

instance Phi Double where
  phi = (sqrt 5 + 1) / 2

------------------------------------------------------------

type R = Float

-- | The 8 vertices of the 16-cell, obtained from (0,0,0,±1) by permuting
-- coordinates.
vs8 :: [Quat R]
vs8 =
  do a <- [1, -1]
     [Q a 0 0 0, Q 0 a 0 0, Q 0 0 a 0, Q 0 0 0 a]

-- | The 16 vertices of the 8-cell, of the form (±½,±½,±½,±½).
vs16 :: [Quat R]
vs16 =
  do w <- [1/2, -1/2]
     x <- [1/2, -1/2]
     y <- [1/2, -1/2]
     z <- [1/2, -1/2]
     return (Q w x y z)

-- | The 24 vertices of the 24-cell, consisting of the union of the
-- vertices of the 8-cell and 16-cell.
vs24 :: [Quat R]
vs24 = vs8 ++ vs16

-- | The 120 vertices of the 600-cell: 16 vertices of the form
-- (±½,±½,±½,±½), and 8 vertices obtained from (0,0,0,±1) by permuting
-- coordinates. The remaining 96 vertices are obtained by taking even
-- permutations of ½(±φ,±1,±1/φ,0).
vs120 :: [Quat R]
vs120 = vs24 ++ vs96
  where
    vs96 =
      do a <- [phi/2, -phi/2]
         b <- [1/2, -1/2]
         c <- [1/phi/2, -1/phi/2]
         x <- [ Q a b c 0
              , Q a c 0 b
              , Q a 0 b c
              , Q b a 0 c --
              , Q b c a 0
              , Q b 0 c a
              , Q c a b 0
              , Q c b 0 a
              , Q c 0 a b --
              , Q 0 a c b
              , Q 0 b a c
              , Q 0 c b a --
              ]
         return x

sizeQ (Q w x y z) = w*w + x*x + y*y + z*z

normalizeQ (Q w x y z) = Q (w/s) (x/s) (y/s) (z/s)
  where s = sqrt (w*w + x*x + y*y + z*z)

-- | The 96 edges of the 24-cell.
edges96 :: [(Quat R, Quat R)]
edges96 = go vs24
  where go [] = []
        go (u : vs) = [ (u, v) | v <- vs, sizeQ (u - v) < 1.25 ] ++ go vs

-- | The 720 edges of the 600-cell.
edges720 :: [(Quat R, Quat R)]
edges720 = go vs120
  where go [] = []
        go (u : vs) = [ (u, v) | v <- vs, sizeQ (u - v) < 0.5 ] ++ go vs
