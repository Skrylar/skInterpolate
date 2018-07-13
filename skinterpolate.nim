
# Cubic interpolation maths are from http://www.paulinternet.nl/?page=bicubic

# XXX i have no doubt that SSE/AVX could be used to speed
# these up quite a lot

proc linear* [T: SomeReal](t, start, stop: T): T {.inline.} =
  ## Bog standard linear interpolation with a time in [0,
  ## 1], between a given `start` and `stop` point.
  return start + (t * (stop - start))

proc linear_delta* [T: SomeReal](t, start, delta: T): T {.inline.} =
  ## Bog standard linear interpolation with a time in
  ## [0, 1], between a given `start` and `start + delta`
  ## point. Slightly more efficient.
  return start + (t * delta)

proc bilinear* [T: SomeReal](x, y, topleft, topright, bottomleft, bottomright: T): T =
  ## Bilinear interpolation between [0, 1], with four known
  ## corners. Right corners are specified as deltas from
  ## the left corner, eliminating one subtraction from some
  ## calculations and being marginally more efficient.
  let top = linear(x, topleft, topright)
  let bottom = linear(x, bottomleft, bottomright)
  return linear(y, top, bottom)

proc bilinear_delta* [T: SomeReal](x, y, topleft, topdelta, bottomleft, bottomdelta: T): T =
  ## Bilinear interpolation between [0, 1], with four known
  ## corners. Right corners are specified as deltas from
  ## the left corner, eliminating one subtraction from some
  ## calculations and being marginally more efficient.
  let top = linear_delta(x, topleft, topdelta)
  let bottom = linear_delta(x, bottomleft, bottomdelta)
  return linear(y, top, bottom)

proc trilinear* [T: SomeReal](x, y, z, fronttopleft, fronttopright, frontbottomleft, frontbottomright, backtopleft, backtopright, backbottomleft, backbottomright: T): T =
  ## Trilinear interpolation between [0, 1], with eight
  ## known corners.
  let fronttop = linear(x, fronttopleft, fronttopright)
  let frontbottom = linear(x, frontbottomleft, frontbottomright)
  let backtop = linear(x, backtopleft, backtopright)
  let backbottom = linear(x, backbottomleft, backbottomright)
  let front = linear(y, fronttop, frontbottom)
  let back = linear(y, backtop, backbottom)
  return linear(z, front, back)

proc trilinear_delta* [T: SomeReal](x, y, z, fronttopleft, fronttopdelta, frontbottomleft, frontbottomdelta, backtopleft, backtopdelta, backbottomleft, backbottomdelta: T): T =
  ## Trilinear interpolation between [0, 1], with eight
  ## known corners. Right corners are specified as deltas
  ## from the left corner, eliminating one subtraction from
  ## some calculations and being marginally more efficient.
  let fronttop = linear_delta(x, fronttopleft, fronttopdelta)
  let frontbottom = linear_delta(x, frontbottomleft, frontbottomdelta)
  let backtop = linear_delta(x, backtopleft, backtopdelta)
  let backbottom = linear_delta(x, backbottomleft, backbottomdelta)
  let front = linear(y, fronttop, frontbottom)
  let back = linear(y, backtop, backbottom)
  return linear(z, front, back)

proc cubic* [T:SomeReal] (t, vm1, v0, v1, v2: T): T =
  ## Performs cubic interpolation over a set of known values.
  ##
  ## - t A: range of [0, 1].
  ## - vm1: The value of f(-1).
  ## - v0: The value of f(0).
  ## - v1: The value of f(1).
  ## - v2: The value of f(2).
  ## 
  ## Cubic interpolation creates a curve over known values,
  ## and uses a derivative to smooth the results. This
  ## derivative is calculated from the values before 0 and
  ## after 1, or "the data point behind the current pair, and
  ## after the current pair." These are stored as vm1 and v2
  ## respectively. If that data is not available (for example,
  ## at the start or end of a list of data points) you may
  ## supply dummy values (ex. supply the value of the first
  ## point as f(-1) and f(0)) or use some interpolated value.

  let a = -((1/2) * vm1) + ((3/2) * v0) - ((3/2) * v1) + ((1/2) * v2)
  let b = vm1 - ((5/2) * v0) + (2 * v1) - ((1/2) * v2)
  let c = (-(1/2) * vm1) + ((1/2) * v1)
  let d = v0

  return (a * t * t * t) +
         (b * t * t) +
         (c * t) +
         d

template catmull_rom* [T:SomeReal] (t, vm1, v0, v1, v2: T): T =
  ## A line created by cubic interpolation, where derivatives
  ## are calculated from additional data points, is known
  ## as a Catmull-Rom spline. It will always pass through
  ## the supplied points, although the rest is smoothed in
  ## some way. Identical to the `cubic` function, but is more
  ## specific about intent.
  cubic[T](t, vm1, v0, v1, v2)

proc cubic* [T:SomeReal] (t, v0, v1: T): T =
  ## Performs cubic interpolation between two values; while
  ## cubic interpolation requires four to do its best job,
  ## the extra points are extrapolated in this call.
  return cubic(t, ((2 * v0) - v1), v0, v1, ((2 * v1) - v0))

proc bicubic* [T:SomeReal](x, y
                           tlm1, tl0, tl1, tl2,
                           trm1, tr0, tr1, tr2,
                           blm1, bl0, bl1, bl2,
                           brm1, br0, br1, br2: T): T =
  ## Performs bicubic interpolation. You must specify the
  ## four necessary points (f(-1), f(0), f(1), and f(20))
  ## for the top left, top right, bottom left and bottom
  ## right of the space.
  let a = cubic(x, tlm1, tl0, tl1, tl2)
  let b = cubic(x, trm1, tr0, tr1, tr2)
  let c = cubic(x, blm1, bl0, bl1, bl2)
  let d = cubic(x, brm1, br0, br1, br2)
  return cubic(y, a, b, c, d)

when isMainModule:
  import unittest, fenv

  suite "Linear interpolation":
    test "Time Zero":
      check abs(linear[float64](0.0, 0, 100) - 0) < epsilon(float64)
    test "Time One":
      check abs(linear[float64](1.0, 0, 100) - 100) < epsilon(float64)

  suite "Cubic interpolation, two points":
    test "Time Zero":
      check abs(cubic[float64](0.0, 0.0, 1.0) - 0.0) < epsilon(float64)
    test "Time One":
      check abs(cubic[float64](1.0, 0.0, 1.0) - 1.0) < epsilon(float64)

  suite "Cubic interpolation, four points":
    test "Time Zero":
      check abs(cubic[float64](0.0, -100.0, 0.0, 1.0, 100.0) - 0.0) < epsilon(float64)
    test "Time One":
      check abs(cubic[float64](1.0, -100.0, 0.0, 1.0, 100.0) - 1.0) < epsilon(float64)

