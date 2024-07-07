{
  /* <script> */
}
/**
 * Single file implementation of sliding windowed infinite Fourier transform (SWIFT)
 *
 * The frequency bands data is formatted like:
 * {lo: lowerBound,
 *  ctr: center,
 *  hi: higherBound}
 *
 * where lo and hi are used for calculating the necessary bandwidth for variable-Q transform spectrum visualizations and ctr for center frequency. This is generated using functions like generateFreqBands
 */
class SWIFT {
  constructor(...args) {
    // initialize the sDFT coefficients
    this.calcCoeffs(args)
    this.spectrumData = []
  }

  calcCoeffs(
    freqBands,
    order = 4,
    timeRes = 600,
    bandwidth = 1,
    sampleRate = 44100,
    compensateBW = true
  ) {
    // calcCoeffs() can be called anywhere else to re-initialize sliding DFT after changes in frequency band distributions and note that x and y are used instead of real and imaginary since vector rotation is the equivalent of the complex one
    this._coeffs = []
    freqBands.map((x, i) => {
      // rX and rY are calculated in advance here since calculating sin and cos functions are pretty slow af
      this._coeffs[i] = {
        rX: Math.cos(((x.ctr * Math.PI) / sampleRate) * 2),
        rY: Math.sin(((x.ctr * Math.PI) / sampleRate) * 2),
        decay:
          Math.E **
          (((-Math.abs(x.hi - x.lo) * Math.PI * bandwidth) / sampleRate -
            1 / ((timeRes * sampleRate) / (Math.PI * 1000))) *
            (compensateBW ? Math.sqrt(order) : 1)),
        coeffs: []
      }
      for (let j = 0; j < order; j++) {
        this._coeffs[i].coeffs[j] = {
          x: 0,
          y: 0
        }
      }
    })
  }

  analyze(dataArray) {
    const newSpectrumData = new Array(this._coeffs.length).fill(0)
    for (const x of dataArray) {
      for (let i = 0; i < this._coeffs.length; i++) {
        for (let j = 0; j < this._coeffs[i].coeffs.length; j++) {
          const input =
              j <= 0
                ? {
                    x: x,
                    y: 0
                  }
                : this._coeffs[i].coeffs[j - 1],
            outX =
              (this._coeffs[i].coeffs[j].x * this._coeffs[i].rX -
                this._coeffs[i].coeffs[j].y * this._coeffs[i].rY) *
                this._coeffs[i].decay +
              input.x * (1 - this._coeffs[i].decay),
            outY =
              (this._coeffs[i].coeffs[j].x * this._coeffs[i].rY +
                this._coeffs[i].coeffs[j].y * this._coeffs[i].rX) *
                this._coeffs[i].decay +
              input.y * (1 - this._coeffs[i].decay)

          this._coeffs[i].coeffs[j].x = outX
          this._coeffs[i].coeffs[j].y = outY
        }
        newSpectrumData[i] = Math.max(
          newSpectrumData[i],
          this._coeffs[i].coeffs[this._coeffs[i].coeffs.length - 1].x ** 2 +
            this._coeffs[i].coeffs[this._coeffs[i].coeffs.length - 1].y ** 2
        )
      }
    }
    this.spectrumData = newSpectrumData.map((x) => Math.sqrt(x))
  }
}
{
  /* </script> */
}
