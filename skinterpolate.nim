
proc linear* [T: SomeReal](t, start, stop: T): T {.inline.} =
  ## Bog standard linear interpolation with a time in [0,
  ## 1], between a given `start` and `stop` point.
  return start + (t * (stop - start))

when isMainModule:
  import unittest, fenv

  suite "Linear interpolation, one dimension"
    test "Time Zero":
      require((linear[float64](0.0, 0, 100) - 0) < epsilon(float64))
    test "Time One":
      require((linear[float64](1.0, 0, 100) - 100) < epsilon(float64))

