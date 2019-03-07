module Quaternion where

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
