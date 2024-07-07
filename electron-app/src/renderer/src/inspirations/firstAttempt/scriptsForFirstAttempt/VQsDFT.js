{
  /* <script> */
}
/**
 * Single file implementation of variable-Q sliding DFT (VQ-sDFT)
 *
 * The frequency bands data is formatted like:
 * {lo: lowerBound,
 *  ctr: center,
 *  hi: higherBound}
 *
 * where lo and hi are used for calculating the necessary bandwidth for variable-Q/constant-Q transform spectrum analysis and ctr for center frequency. This is generated using functions like generateFreqBands()
 *
 * Note: This algorithm is derived from the paper "Application of Improved Sliding DFT Algorithm for Non-Integer k" by Carl Q. Howard (https://acoustics.asn.au/conference_proceedings/AAS2021/papers/p60.pdf)
 */
export class VQsDFT {
  constructor(...args) {
    this.calcCoeffs(args)
    this.spectrumData = []
  }

  calcCoeffs(
    freqBands,
    window = [1, 0.5],
    timeRes = 600,
    bandwidth = 1,
    bufferSize = 44100,
    sampleRate = 44100,
    useNC = false
  ) {
    this._coeffs = freqBands.map((x) => {
      const fiddles = [],
        twiddles = [],
        resonCoeffs = [],
        coeffs1 = [],
        coeffs2 = [],
        coeffs3 = [],
        coeffs4 = [],
        coeffs5 = [],
        gains = [],
        period = Math.trunc(
          Math.min(
            bufferSize,
            sampleRate / (bandwidth * Math.abs(x.hi - x.lo) + 1 / (timeRes / 1000))
          )
        ), // N must be an integer, but K doesn't have to be
        minIdx = useNC ? 0 : -window.length + 1,
        maxIdx = useNC ? 2 : window.length
      // this below is needed since we have to apply a frequency-domain window function
      for (let i = minIdx; i < maxIdx; i++) {
        const amplitude = useNC ? 1 : window[Math.abs(i)] * (-(Math.abs(i) % 2) * 2 + 1),
          k = (x.ctr * period) / sampleRate + i - useNC / 2,
          fid = -2 * Math.PI * k,
          twid = (2 * Math.PI * k) / period,
          reson = 2 * Math.cos((2 * Math.PI * k) / period)
        fiddles.push({
          x: Math.cos(fid),
          y: Math.sin(fid)
        })
        twiddles.push({
          x: Math.cos(twid),
          y: Math.sin(twid)
        })
        resonCoeffs.push(reson)
        coeffs1.push({ x: 0, y: 0 })
        coeffs2.push({ x: 0, y: 0 })
        coeffs3.push({ x: 0, y: 0 })
        coeffs4.push({ x: 0, y: 0 })
        coeffs5.push({ x: 0, y: 0 })
        gains.push(amplitude)
      }
      return {
        period: period,
        twiddles: twiddles,
        fiddles: fiddles,
        resonCoeffs: resonCoeffs,
        coeffs1: coeffs1,
        coeffs2: coeffs2,
        coeffs3: coeffs3,
        coeffs4: coeffs4,
        coeffs5: coeffs5,
        gains: gains,
        nc: useNC
      }
    })
    this._buffer = new Array(bufferSize + 1).fill(0)
    this._bufferIdx = this._buffer.length - 1 // this is required for circular buffer
  }

  analyze(samples) {
    this.spectrumData = new Array(this._coeffs.length).fill(0)
    for (const sample of samples) {
      // Admittedly slow linear buffer
      /*
      this._buffer.push(sample);
      this._buffer.shift();
      */
      // Circular buffer
      this._bufferIdx =
        (((this._bufferIdx + 1) % this._buffer.length) + this._buffer.length) % this._buffer.length
      this._buffer[this._bufferIdx] = sample
      for (let i = 0; i < this._coeffs.length; i++) {
        const coeff = this._coeffs[i],
          kernelLength = coeff.coeffs1.length,
          /*oldest = this._buffer.length-coeff.period-1,
              latest = this._buffer.length-1,*/
          oldest =
            (((this._bufferIdx - coeff.period) % this._buffer.length) + this._buffer.length) %
            this._buffer.length,
          latest = this._bufferIdx,
          sum = {
            x: 0,
            y: 0
          }
        for (let j = 0; j < kernelLength; j++) {
          const fiddle = coeff.fiddles[j],
            twiddle = coeff.twiddles[j],
            // Comb stage
            combX = this._buffer[latest] * fiddle.x - this._buffer[oldest],
            combY = this._buffer[latest] * fiddle.y

          // Second stage
          coeff.coeffs1[j].x = combX * twiddle.x - combY * twiddle.y - coeff.coeffs2[j].x
          coeff.coeffs1[j].y = combX * twiddle.y + combY * twiddle.x - coeff.coeffs2[j].y

          coeff.coeffs2[j].x = combX
          coeff.coeffs2[j].y = combY

          // Real resonator
          coeff.coeffs3[j].x =
            coeff.coeffs1[j].x + coeff.resonCoeffs[j] * coeff.coeffs4[j].x - coeff.coeffs5[j].x
          coeff.coeffs3[j].y =
            coeff.coeffs1[j].y + coeff.resonCoeffs[j] * coeff.coeffs4[j].y - coeff.coeffs5[j].y

          coeff.coeffs5[j].x = coeff.coeffs4[j].x
          coeff.coeffs5[j].y = coeff.coeffs4[j].y

          coeff.coeffs4[j].x = coeff.coeffs3[j].x
          coeff.coeffs4[j].y = coeff.coeffs3[j].y

          sum.x += (coeff.coeffs3[j].x * coeff.gains[j]) / coeff.period
          sum.y += (coeff.coeffs3[j].y * coeff.gains[j]) / coeff.period
        }
        const period = coeff.period
        this.spectrumData[i] = Math.max(
          this.spectrumData[i],
          coeff.nc
            ? -(((coeff.coeffs3[0].x / period) * coeff.coeffs3[1].x) / period) -
                ((coeff.coeffs3[0].y / period) * coeff.coeffs3[1].y) / period
            : sum.x ** 2 + sum.y ** 2
        )
      }
    }
    this.spectrumData = this.spectrumData.map((x) => Math.sqrt(x))
  }
}
{
  /* </script> */
}
