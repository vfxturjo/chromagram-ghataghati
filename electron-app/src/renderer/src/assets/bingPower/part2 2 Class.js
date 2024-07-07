// chromagram.js
import FFT from 'fftw-js/src/main'

class ChromagramCalculator {
  constructor() {
    this.NUM_BINS = 12
    this.SAMPLE_RATE = 44100
    this.FFT_SIZE = 4096
    this.FFTWorker = FFT(this.FFT_SIZE)
    this.BUFFER_LEN = this.FFT_SIZE / 2
    this.circularBuffer = new Float32Array(this.BUFFER_LEN)
    this.bufferIndex = 0
    this.selectedWindow = 'Hanning' // Default window function
    this.sigma = 0.4 // Default sigma for Gaussian window
  }

  addToBuffer(samples) {
    samples.forEach((sample) => {
      this.circularBuffer[this.bufferIndex] = sample
      this.bufferIndex = (this.bufferIndex + 1) % this.BUFFER_LEN
    })
  }

  calculateChromagram() {
    const windowedSamples = new Float32Array(this.FFT_SIZE)

    for (let i = 0; i < this.FFT_SIZE; i++) {
      const circularIndex = (this.bufferIndex - i + this.BUFFER_LEN) % this.BUFFER_LEN
      windowedSamples[i] = this.circularBuffer[circularIndex] * this.windowFunction(i)
    }

    const frequencyData = this.FFTWorker.forward(windowedSamples)

    const chromagramValues = new Array(this.NUM_BINS).fill(0)
    frequencyData.forEach((magnitude, bin) => {
      if (bin < this.FFT_SIZE / 2) {
        const frequency = (bin * this.SAMPLE_RATE) / this.FFT_SIZE
        const pitchClass = Math.round(12 * Math.log2(frequency / 440)) % 12
        chromagramValues[pitchClass] += magnitude
      }
    })

    const totalEnergy = chromagramValues.reduce((sum, val) => sum + val, 0)
    if (totalEnergy > 0) {
      chromagramValues.forEach((val, idx) => {
        chromagramValues[idx] = val / totalEnergy
      })
    }

    return chromagramValues
  }

  windowFunction(index) {
    if (this.selectedWindow === 'Hanning') {
      return 0.5 * (1 - Math.cos((2 * Math.PI * index) / (this.FFT_SIZE - 1)))
    } else if (this.selectedWindow === 'Gaussian') {
      const center = (this.FFT_SIZE - 1) / 2
      return Math.exp(-0.5 * ((index - center) / (this.sigma * center)) ** 2)
    } else {
      // Default to rectangular window
      return 1
    }
  }

  setWindowFunction(windowName) {
    this.selectedWindow = windowName
  }

  setSigma(sigma) {
    this.sigma = sigma
  }
}

export default ChromagramCalculator
