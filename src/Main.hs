module Main where

import Graphics.Gloss
import Graphics.Gloss.Interface.Pure.Game

{-
------------------------------------------------------------
Notes on the geometry of a 3-sphere.
------------------------------------------------------------

Points in space are described by 4-tuples (w,x,y,z) in R^4 such that
w^2 + x^2 + y^2 + z^2 = 1.

We take point o = (1,0,0,0) as the origin, or starting point.

At each point p on the 3-sphere, we assign vectors E(p), N(p), and
U(p) that identify the directions East, North, and Up. All three
vectors must be orthogonal to each other, and to p, for any p. All
three mappings should also be linear.

At the origin, we have:
E(1,0,0,0) = (0,1,0,0)
N(1,0,0,0) = (0,0,1,0)
U(1,0,0,0) = (0,0,0,1)

Right-handed direction mapping:
I(w,x,y,z) = (+w, +x, +y, +z)
E(w,x,y,z) = (-x, +w, +z, -y)
N(w,x,y,z) = (-y, -z, +w, +x)
U(w,x,y,z) = (-z, +y, -x, +w)

Unitary, direction-preserving transformations:

Origin (1,0,0,0) to east pole (0,1,0,0).
Just east of origin  (1,x,0,0) to just east of pole  (-x,1,0,0).
Just north of origin (1,0,y,0) to just north of pole (0,1,0,y).
Just above origin    (1,0,0,z) to just above pole    (0,1,-z,0).

oe(w,x,y,z) = (-x,+w,-z,+y)

[ 0 -1  0  0] [w]   [-x]
[+1  0  0  0] [x] = [+w]
[ 0  0  0 -1] [y]   [-z]
[ 0  0 +1  0] [z]   [+y]


Origin (1,0,0,0) to point p = (a,b,c,d).
Just east of origin  (1,x,0,0)
to just east of p    (a,b,c,d) + x*E(a,b,c,d).
                     (a-bx, b+ax, c+dx, d-cx)

Just north of origin (1,0,y,0)
to just north of p   (a,b,c,d) + y*N(a,b,c,d).
                     (a-cy, b-dy, c+ay, d+by)

Just above origin    (1,0,0,z)
to just above p      (a,b,c,d) + z*U(a,b,c,d.
                     (a-dz, b+cz, c-bz, d+az)


T{a,b,c,d}(w,x,y,z) = (aw-bx-cy-dz, bw+ax-dy+cz, cw+dx+ay-bz, dw-cx+by+az)


[+a -b -c -d] [w]   [aw-bx-cy-dz]
[+b +a -d +c] [x] = [bw+ax-dy+cz]
[+c +d +a -b] [y]   [cw+dx+ay-bz]
[+d -c +b +a] [z]   [dw-cx+by+az]

Determinant of T is (a^2 + b^2 + c^2 + d^2)^2 = 1^2 = 1.

    [+a -b -c -d]
T = [+b +a -d +c]
    [+c +d +a -b]
    [+d -c +b +a]

    [ 0 -1  0  0]
E = [+1  0  0  0]
    [ 0  0  0 +1]
    [ 0  0 -1  0]

    [ 0  0 -1  0]
N = [ 0  0  0 -1]
    [+1  0  0  0]
    [ 0 +1  0  0]

    [ 0  0  0 -1]
U = [ 0  0 +1  0]
    [ 0 -1  0  0]
    [+1  0  0  0]

Each commutes with T, i.e. TE = ET, etc.

T*transp(T) =

[+a -b -c -d]   [+a +b +c +d]   [1 0 0 0]
[+b +a -d +c] * [-b +a +d -c] = [0 1 0 0]
[+c +d +a -b]   [-c -d +a +b]   [0 0 1 0]
[+d -c +b +a]   [-d +c -b +a]   [0 0 0 1]

------------------------------------------------------------

The 120 vertices of the 600-cell:

16 vertices of the form (±½,±½,±½,±½), and 8 vertices obtained from
(0,0,0,±1) by permuting coordinates. The remaining 96 vertices are
obtained by taking even permutations of ½(±φ,±1,±1/φ,0).
-}

import Data.List ((\\), sortBy)

data Quat a = Q a a a a
  deriving (Eq, Ord, Show)

wQ, xQ, yQ, zQ :: Quat a -> a
wQ (Q w _ _ _) = w
xQ (Q _ x _ _) = x
yQ (Q _ _ y _) = y
zQ (Q _ _ _ z) = z

instance Num a => Num (Quat a) where
  Q a b c d + Q w x y z = Q (a+w) (b+x) (c+y) (d+z)
  Q a b c d - Q w x y z = Q (a-w) (b-x) (c-y) (d-z)
  Q a b c d * Q w x y z = Q (a*w - b*x - c*y - d*z)
                            (b*w + a*x - d*y + c*z)
                            (c*w + d*x + a*y - b*z)
                            (d*w - c*x + b*y + a*z)
  negate (Q a b c d) = Q (negate a) (negate b) (negate c) (negate d)
  fromInteger n = Q (fromInteger n) 0 0 0
  abs _ = error "abs Quat"
  signum _ = error "signum Quat"

conjQ :: Num a => Quat a -> Quat a
conjQ (Q a b c d) = Q a (-b) (-c) (-d)

instance Fractional a => Fractional (Quat a) where
  recip (Q a b c d) = let s = a*a + b*b + c*c + d*d
                      in Q (a/s) (-b/s) (-c/s) (-d/s)
  fromRational r = Q (fromRational r) 0 0 0

class Phi a where
  phi :: a

instance Phi Float where
  phi = (sqrt 5 + 1) / 2

instance Phi Double where
  phi = (sqrt 5 + 1) / 2

-- | The idea is that turn(a,b,c,d) should rotate (1,0,0,0) onto
-- (a,b,c,d), while preserving directions for all other points. It
-- turns out that this is precisely quaternion multiplication, where
-- (a,b,c,d) represents a + bi + cj + dk.

turn (a,b,c,d) (w,x,y,z) =
  (a*w - b*x - c*y - d*z,
   b*w + a*x - d*y + c*z,
   c*w + d*x + a*y - b*z,
   d*w - c*x + b*y + a*z)

unturn (a,b,c,d) (w,x,y,z) =
  (a*w + b*x + c*y + d*z,
  -b*w + a*x + d*y - c*z,
  -c*w - d*x + a*y + b*z,
  -d*w + c*x - b*y + a*z)

data P = P Rational Rational
  deriving (Eq, Show)

instance Ord P where
  compare (P a b) (P c d) =
      case compare a c of
        EQ             -> compare b d
        LT | b <= d    -> LT
           | otherwise -> let x = (b-d)/(c-a)
                          in compare (x*x) (x+1)
        GT | b >= d    -> GT
           | otherwise -> let x = (d-b)/(a-c)
                          in compare (x+1) (x*x)

instance Num P where
  P a b + P c d = P (a + c) (b + d)
  P a b - P c d = P (a - c) (b - d)
  P a b * P c d = P (a*c + b*c + a*d) (a*c + b*d)
  negate (P a b) = P (negate a) (negate b)
  fromInteger n = P 0 (fromInteger n)
  abs (P _ _) = error "abs P"
  signum (P _ _) = error "signum P"

instance Fractional P where
  recip (P a b) = P (-a/d) ((a+b)/d)
      where d = b*b + a*b - a*a
  fromRational r = P 0 r

instance Phi P where
  phi = P 1 0

unP (P a b) = fromRational a * phi + fromRational b

{-
16 vertices of the form (±½,±½,±½,±½), and 8 vertices obtained from
(0,0,0,±1) by permuting coordinates. The remaining 96 vertices are
obtained by taking even permutations of ½(±φ,±1,±1/φ,0).
-}

type R = Float

vs8 :: [Quat R]
vs8 = do
  a <- [1, -1]
  [Q a 0 0 0, Q 0 a 0 0, Q 0 0 a 0, Q 0 0 0 a]

vs24 :: [Quat R]
vs24 = vs8 ++ vs16
  where
    vs16 = do
      w <- [1/2, -1/2]
      x <- [1/2, -1/2]
      y <- [1/2, -1/2]
      z <- [1/2, -1/2]
      return (Q w x y z)

vs120 :: [Quat R]
vs120 = vs24 ++ vs96
  where
    vs96 = do
      a <- [phi/2, -phi/2]
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

edges720 :: [(Quat R, Quat R)]
edges720 = go vs120
  where go [] = []
        go (u : vs) = [ (u, v) | v <- vs, sizeQ (u - v) < 0.5 ] ++ go vs

edges96 :: [(Quat R, Quat R)]
edges96 = go vs24
  where go [] = []
        go (u : vs) = [ (u, v) | v <- vs, sizeQ (u - v) < 1.25 ] ++ go vs

edgeNodes :: Int -> (Quat R, Quat R) -> [Quat R]
edgeNodes n (u, v) = map normalizeQ [ u * d + v * (1 - d) | d <- ds ]
  where
    ds = [ fromIntegral a / fromIntegral n | a <- [1 .. n-1] ]

someEdges :: [(Quat R, Quat R)]
someEdges = [ (t * u1, t * u2) | t <- vs120, (u1, u2) <- es ]
    where
      es = [(v1, v2), (v2, v3), (v3, v1)]
      v1 = vs120 !! 0
      v2 = vs120 !! 38
      v3 = vs120 !! 48


allNodes = vs120 ++ concatMap (edgeNodes 20) edges720
--allNodes = vs24 ++ concatMap (edgeNodes 20) edges

type Ball = (Quat R, Color, Float, Int)

allBalls :: [Ball]
allBalls =
  [ (v, blue, 0.08, 0) | v <- take 1 vs120 ] ++
  [ (v, red, 0.08, 0) | v <- drop 1 $ take 24 vs120 ] ++
  [ (v, black, 0.025, 0) | v <- drop 24 $ vs120 ] ++
  [ (v, black, 0.02, 0) | v <- concatMap (edgeNodes 24) edges720 ] ++
  --[ (v, black, 0.02, 0) | v <- concatMap (edgeNodes 24) someEdges ] ++
  --[ (v, black, 0.025, 0) | v <- concatMap (edgeNodes 32) edges96 ] ++
  --[ (t * u, c, s, n1*n2) | (t, n1) <- zip vs24 [1..], (u, c, s, n2) <- flag ]
  [ (t * u, c, s, n1*n2) | (t, n1) <- zip vs120 [1..], (u, c, s, n2) <- flag ]

flag :: [Ball]
flag =
  [ (normalizeQ $ Q 1 0 0 0.04, red, 0.025, 0)
  , (normalizeQ $ Q 1 0 0 0.08, red, 0.025, 0)
  , (normalizeQ $ Q 1 0 0 0.12, red, 0.025, 0)
  , (normalizeQ $ Q 1 0 0 0.16, red, 0.025, 0)
  , (normalizeQ $ Q 1 0.04 0 0.14, blue, 0.025, 1)
  ]



{-
play
:: forall world .
=> Display    Display mode.
-> Color      Background color.
-> Int	      Number of simulation steps to take for each second of real time.
-> world      The initial world.
-> (world -> Picture)	A function to convert the world a picture.
-> (Event -> world -> world)	A function to handle input events.
-> (Float -> world -> world)	A function to step the world one iteration.
        It is passed the period of time (in seconds) needing to be advanced.
-> IO ()	 Play a game in a window.
-}

addR3 :: (R, R, R) -> (R, R, R) -> (R, R, R)
addR3 (x1, y1, z1) (x2, y2, z2) = (x1 + x2, y1 + y2, z1 + z2)

subR3 :: (R, R, R) -> (R, R, R) -> (R, R, R)
subR3 (x1, y1, z1) (x2, y2, z2) = (x1 - x2, y1 - y2, z1 - z2)

data World = World { position :: Quat R, velocity :: (R, R, R) }

worldInit :: World
worldInit = World { position = Q 1 0 0 0, velocity = (0, 0, 0) }

advanceWorld :: Float -> World -> World
advanceWorld t world = world { position = normalizeQ (p * Q 1 dx dy dz) }
  where
    p = position world
    (vx, vy, vz) = velocity world
    (dx, dy, dz) = (t*vx, t*vy, t*vz)

eventWorld :: Event -> World -> World
eventWorld (EventKey (Char c) ks _ _) world =
  case ks of
    Down -> world { velocity = velocity world `addR3` vel c }
    Up   -> world { velocity = velocity world `subR3` vel c }
  where
    p = 0.5 -- maximum speed
    m = -p
    vel 'i' = (0, p, 0)
    vel 'j' = (m, 0, 0)
    vel 'k' = (0, m, 0)
    vel 'l' = (p, 0, 0)
    vel 'u' = (0, 0, p)
    vel 'o' = (0, 0, m)
    vel _   = (0, 0, 0)
eventWorld _ world = world


filledCircle :: Float -> Picture
filledCircle x = ThickCircle (x/2) x



drawWorld :: World -> Picture
drawWorld world =
  Pictures
  [ Scale 200 200 $ Pictures (map showBall (reverse (sortBy comp ps)))
  , Translate (-300) (-200) $ Scale 20 4 $ Pictures $
    [ Color red   $ Polygon [(0,0),(0,1),(a,1),(a,0)]
    , Color blue  $ Polygon [(0,2),(0,3),(b,3),(b,2)]
    , Color green $ Polygon [(0,4),(0,5),(c,5),(c,4)]
    , Color black $ Polygon [(0,6),(0,7),(d,7),(d,6)]
    ]
  , Translate (-300) (-250) $ Scale (0.1) (0.1) $ Text (show (length ps))
  , Line [(-4.5, 0), (4.5, 0)]
  , Line [(0, -4.5), (0, 4.5)]
  ]
  where
    Q a b c d = position world
    comp (a, _, _, _) (b, _, _, _) = compare (yQ a) (yQ b)
    t = conjQ (position world)
    vs' = [ (t * v, c, s, n) | (v, c, s, n) <- allBalls ]
    ps = filter (\(Q w x y z, _, _, _) -> w > 0.01 && y > 0.02) vs'
    showBall (Q w x y z, c, s, n) =
--      Translate x z $
--      filledCircle 0.05
      Translate (x/y) (z/y) $
      Pictures [ Color (gradientColor (w*w) c bgColor) $ filledCircle (s / y)
               , Color white $
                   Scale (s / y / 200) (s / y / 200) $
                   if n > 0 then Text (show n) else Blank
               ]

bgColor :: Color
bgColor = azure `addColors` white `addColors` white

gradientColor x = mixColors x (1 - x)

main :: IO ()
main = play
  (InWindow "Three-Sphere" (800, 600) (10, 10))
  bgColor
  30 --fps
  worldInit
  drawWorld
  eventWorld
  advanceWorld


{-
   | E  W  N  S  U  D
---+-----------------
+w |+x -x +y -y +z -z
-w |
+x |-w +w +z -z -y +y
-x |
+y |-z +z -w +w +x -x
-y |
+z |+y -y -x +x -w +w
-z |

E/W:
+w +x -w -x
+z +y -z -y

N/S:
+w +y -w -y
+x +z -x -z

U/D:
+w +z -w -z
+y +x -y -x

-}
