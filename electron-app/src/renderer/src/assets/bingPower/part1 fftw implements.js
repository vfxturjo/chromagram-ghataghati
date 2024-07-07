// chromagram.js
import FFT from 'fftw-js/src/main'

const NUM_BINS = 12
const SAMPLE_RATE = 44100
const FFT_SIZE = 2048
const BUFFER_LEN = FFT_SIZE / 2

const circularBuffer = new Float32Array(BUFFER_LEN)
let bufferIndex = 0

function addToBuffer(samples) {
  samples.forEach((sample) => {
    circularBuffer[bufferIndex] = sample
    bufferIndex = (bufferIndex + 1) % BUFFER_LEN
  })
}

function calculateChromagram() {
  const windowedSamples = new Float32Array(FFT_SIZE)

  for (let i = 0; i < FFT_SIZE; i++) {
    const circularIndex = (bufferIndex + i) % BUFFER_LEN
    windowedSamples[i] = circularBuffer[circularIndex]
  }

  // Use FFTW-JS for FFT
  const fftResult = FFT.FFT(windowedSamples)

  // Extract magnitude data (you may need to adjust this based on the library)
  const frequencyData = new Float32Array(FFT_SIZE)
  for (let i = 0; i < FFT_SIZE; i++) {
    frequencyData[i] = Math.sqrt(fftResult[i][0] ** 2 + fftResult[i][1] ** 2)
  }

  const chromagramValues = new Array(NUM_BINS).fill(0)
  frequencyData.forEach((magnitude, bin) => {
    const frequency = (bin * SAMPLE_RATE) / FFT_SIZE
    const pitchClassIndex = Math.floor((12 * Math.log2(frequency / 440)) % 12)
    chromagramValues[pitchClassIndex] += magnitude
  })

  // Normalize the chromagram (optional)
  const totalEnergy = chromagramValues.reduce((sum, val) => sum + val, 0)
  chromagramValues.forEach((val, idx) => {
    chromagramValues[idx] = val / totalEnergy
  })

  return chromagramValues
}

export { addToBuffer, calculateChromagram }
