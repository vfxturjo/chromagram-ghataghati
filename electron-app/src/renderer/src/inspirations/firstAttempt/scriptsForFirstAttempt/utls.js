{
  /* <script> */
}

export {
  map,
  clamp,
  idxWrapOver,
  hertzToFFTBin,
  fftBinToHertz,
  calcFFT,
  calcComplexFFT,
  calcComplexInputFFT,
  convolveReal
}

function map(x, min, max, targetMin, targetMax) {
  return ((x - min) / (max - min)) * (targetMax - targetMin) + targetMin
}

function clamp(x, min, max) {
  return Math.min(Math.max(x, min), max)
}

function idxWrapOver(x, length) {
  return ((x % length) + length) % length
}
// Hz and FFT bin conversion
function hertzToFFTBin(x, y = 'round', bufferSize = 4096, sampleRate = 44100) {
  const bin = (x * bufferSize) / sampleRate
  let func = y

  if (!['floor', 'ceil', 'trunc'].includes(func)) func = 'round' // always use round if you specify an invalid/undefined value

  return Math[func](bin)
}

function fftBinToHertz(x, bufferSize = 4096, sampleRate = 44100) {
  return (x * sampleRate) / bufferSize
}

// Calculate the FFT
function calcFFT(input) {
  let fft = input.map((x) => x)
  let fft2 = input.map((x) => x)
  transform(fft, fft2)
  let output = new Array(Math.round(fft.length / 2)).fill(0)
  for (let i = 0; i < output.length; i++) {
    output[i] = Math.hypot(fft[i], fft2[i]) / fft.length
  }
  return output
}

function calcComplexFFT(input) {
  let fft = input.map((x) => x)
  let fft2 = input.map((x) => x)
  transform(fft, fft2)
  return input.map((_, i, arr) => {
    return {
      re: fft[i] / (arr.length / 2),
      im: fft2[i] / (arr.length / 2),
      magnitude: Math.hypot(fft[i], fft2[i]) / (arr.length / 2),
      phase: Math.atan2(fft2[i], fft[i])
    }
  })
}

function calcComplexInputFFT(real, imag) {
  if (real.length !== imag.length) return []
  const fft1 = real.map((x) => x),
    fft2 = imag.map((x) => x)
  transform(fft1, fft2)
  return real.map((_, i, arr) => {
    return {
      re: fft1[i] / arr.length,
      im: fft2[i] / arr.length,
      magnitude: Math.hypot(fft1[i], fft2[i]) / arr.length,
      phase: Math.atan2(fft2[i], fft1[i])
    }
  })
}

/**
 * FFT and convolution (JavaScript)
 *
 * Copyright (c) 2017 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 */

/*
 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
 * The vector can have any length. This is a wrapper function.
 */
function transform(real, imag) {
  const n = real.length
  if (n != imag.length) throw 'Mismatched lengths'
  if (n <= 0) return
  else if (2 ** Math.trunc(Math.log2(n)) === n)
    // Is power of 2
    transformRadix2(real, imag)
  // More complicated algorithm for arbitrary sizes
  else transformBluestein(real, imag)
}

/*
 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
 */
function inverseTransform(real, imag) {
  transform(imag, real)
}

/*
 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
 */
function transformRadix2(real, imag) {
  // Length variables
  const n = real.length
  if (n != imag.length) throw 'Mismatched lengths'
  if (n <= 1)
    // Trivial transform
    return
  const logN = Math.log2(n)
  if (2 ** Math.trunc(logN) !== n) throw 'Length is not a power of 2'

  // Trigonometric tables
  let cosTable = new Array(n / 2)
  let sinTable = new Array(n / 2)
  for (let i = 0; i < n / 2; i++) {
    cosTable[i] = Math.cos((2 * Math.PI * i) / n)
    sinTable[i] = Math.sin((2 * Math.PI * i) / n)
  }

  // Bit-reversed addressing permutation
  for (let i = 0; i < n; i++) {
    let j = reverseBits(i, logN)
    if (j > i) {
      let temp = real[i]
      real[i] = real[j]
      real[j] = temp
      temp = imag[i]
      imag[i] = imag[j]
      imag[j] = temp
    }
  }

  // Cooley-Tukey decimation-in-time radix-2 FFT
  for (let size = 2; size <= n; size *= 2) {
    let halfsize = size / 2
    let tablestep = n / size
    for (let i = 0; i < n; i += size) {
      for (let j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
        const l = j + halfsize
        const tpre = real[l] * cosTable[k] + imag[l] * sinTable[k]
        const tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k]
        real[l] = real[j] - tpre
        imag[l] = imag[j] - tpim
        real[j] += tpre
        imag[j] += tpim
      }
    }
  }

  // Returns the integer whose value is the reverse of the lowest 'bits' bits of the integer 'x'.
  function reverseBits(x, bits) {
    let y = 0
    for (let i = 0; i < bits; i++) {
      y = (y << 1) | (x & 1)
      x >>>= 1
    }
    return y
  }
}

/*
 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
 * Uses Bluestein's chirp z-transform algorithm.
 */
function transformBluestein(real, imag) {
  // Find a power-of-2 convolution length m such that m >= n * 2 + 1
  const n = real.length
  if (n != imag.length) throw 'Mismatched lengths'
  const m = 2 ** Math.trunc(Math.log2(n * 2) + 1)

  // Trignometric tables
  let cosTable = new Array(n)
  let sinTable = new Array(n)
  for (let i = 0; i < n; i++) {
    let j = (i * i) % (n * 2) // This is more accurate than j = i * i
    cosTable[i] = Math.cos((Math.PI * j) / n)
    sinTable[i] = Math.sin((Math.PI * j) / n)
  }

  // Temporary vectors and preprocessing
  let areal = newArrayOfZeros(m)
  let aimag = newArrayOfZeros(m)
  for (let i = 0; i < n; i++) {
    areal[i] = real[i] * cosTable[i] + imag[i] * sinTable[i]
    aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i]
  }
  let breal = newArrayOfZeros(m)
  let bimag = newArrayOfZeros(m)
  breal[0] = cosTable[0]
  bimag[0] = sinTable[0]
  for (let i = 1; i < n; i++) {
    breal[i] = breal[m - i] = cosTable[i]
    bimag[i] = bimag[m - i] = sinTable[i]
  }

  // Convolution
  let creal = new Array(m)
  let cimag = new Array(m)
  convolveComplex(areal, aimag, breal, bimag, creal, cimag)

  // Postprocessing
  for (let i = 0; i < n; i++) {
    real[i] = creal[i] * cosTable[i] + cimag[i] * sinTable[i]
    imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i]
  }
}

/*
 * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
 */
function convolveReal(x, y, out) {
  const n = x.length
  if (n != y.length || n != out.length) throw 'Mismatched lengths'
  convolveComplex(x, newArrayOfZeros(n), y, newArrayOfZeros(n), out, newArrayOfZeros(n))
}

/*
 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
 */
function convolveComplex(xreal, ximag, yreal, yimag, outreal, outimag) {
  const n = xreal.length
  if (
    n != ximag.length ||
    n != yreal.length ||
    n != yimag.length ||
    n != outreal.length ||
    n != outimag.length
  )
    throw 'Mismatched lengths'

  xreal = xreal.slice()
  ximag = ximag.slice()
  yreal = yreal.slice()
  yimag = yimag.slice()
  transform(xreal, ximag)
  transform(yreal, yimag)

  for (let i = 0; i < n; i++) {
    const temp = xreal[i] * yreal[i] - ximag[i] * yimag[i]
    ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i]
    xreal[i] = temp
  }
  inverseTransform(xreal, ximag)

  for (let i = 0; i < n; i++) {
    // Scaling (because this FFT implementation omits it)
    outreal[i] = xreal[i] / n
    outimag[i] = ximag[i] / n
  }
}

function newArrayOfZeros(n) {
  let result = new Array(n).fill(0)
  return result
}
// </script>
