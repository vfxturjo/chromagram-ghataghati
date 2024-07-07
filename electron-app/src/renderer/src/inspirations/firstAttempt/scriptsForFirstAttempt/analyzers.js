{
  /* <script> */
}
class AnalogStyleAnalyzer {
  constructor(...args) {
    // initialize the sDFT coefficients
    this.calcCoeffs(args)
    this.spectrumData = []
  }

  calcCoeffs(
    freqBands,
    order = 4,
    timeRes = Infinity,
    bandwidth = 1,
    sampleRate = 44100,
    compensateBW = true,
    prewarpQ = false
  ) {
    this._coeffs = freqBands.map((x) => {
      // biquad bandpass filter (cascaded biquad bandpass is not Butterworth nor Bessel, rather it is something called "critically-damped" since each filter stage shares the same every biquad coefficients)
      const K = Math.tan((Math.PI * x.ctr) / sampleRate),
        bw = Math.abs(x.hi - x.lo) * bandwidth + 1 / (timeRes / 1000),
        qCompensationFactor = prewarpQ ? (Math.PI * x.ctr) / sampleRate / K : 1,
        Q = ((x.ctr / bw) * qCompensationFactor) / (compensateBW ? Math.sqrt(order) : 1),
        norm = 1 / (1 + K / Q + K * K),
        a0 = (K / Q) * norm,
        a1 = 0,
        a2 = -a0,
        b1 = 2 * (K * K - 1) * norm,
        b2 = (1 - K / Q + K * K) * norm,
        zs = []
      for (let i = 0; i < order; i++) {
        zs[i] = {
          z1: 0,
          z2: 0,
          out: 0
        }
      }
      return {
        a0: a0,
        a1: a1,
        a2: a2,
        b1: b1,
        b2: b2,
        zs: zs
      }
    })
  }

  analyze(samples) {
    const newSpectrumData = new Array(this._coeffs.length).fill(0)
    for (const x of samples) {
      for (let i = 0; i < this._coeffs.length; i++) {
        for (let j = 0; j < this._coeffs[i].zs.length; j++) {
          const input = j <= 0 ? x : this._coeffs[i].zs[j - 1].out
          this._coeffs[i].zs[j].out = input * this._coeffs[i].a0 + this._coeffs[i].zs[j].z1
          this._coeffs[i].zs[j].z1 =
            input * this._coeffs[i].a1 +
            this._coeffs[i].zs[j].z2 -
            this._coeffs[i].b1 * this._coeffs[i].zs[j].out
          this._coeffs[i].zs[j].z2 =
            input * this._coeffs[i].a2 - this._coeffs[i].b2 * this._coeffs[i].zs[j].out
        }
        newSpectrumData[i] = Math.max(
          newSpectrumData[i],
          Math.abs(this._coeffs[i].zs[this._coeffs[i].zs.length - 1].out)
        )
      }
    }
    this.spectrumData = newSpectrumData.map((x) => x / 2)
  }
}
{
  /* </script> */
}
