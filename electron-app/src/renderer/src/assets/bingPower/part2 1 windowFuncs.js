// chromagram.js
import FFT from 'fftw-js/src/main'

// Constants for chromagram calculation
const NUM_BINS = 12
const SAMPLE_RATE = 44100
const FFT_SIZE = 4096
const FFTWorker = FFT(FFT_SIZE)
const BUFFER_LEN = FFT_SIZE / 2

// Circular buffer to store recent audio samples
const circularBuffer = new Float32Array(BUFFER_LEN)
let bufferIndex = 0

// Window function options
const WINDOW_FUNCTIONS = {
  HANNING: 'Hanning',
  GAUSSIAN: 'Gaussian'
}

let selectedWindow = WINDOW_FUNCTIONS.HANNING // Default window function

// Function to add samples to the circular buffer
function addToBuffer(samples) {
  samples.forEach((sample) => {
    circularBuffer[bufferIndex] = sample
    bufferIndex = (bufferIndex + 1) % BUFFER_LEN
  })
}

// Function to calculate the chromagram from the circular buffer
function calculateChromagram() {
  const windowedSamples = new Float32Array(FFT_SIZE)

  for (let i = 0; i < FFT_SIZE; i++) {
    const circularIndex = (bufferIndex - i + BUFFER_LEN) % BUFFER_LEN
    windowedSamples[i] = circularBuffer[circularIndex] * windowFunction(i)
  }

  const frequencyData = FFTWorker.forward(windowedSamples)

  const chromagramValues = new Array(NUM_BINS).fill(0)
  frequencyData.forEach((magnitude, bin) => {
    if (bin < FFT_SIZE / 2) {
      const frequency = (bin * SAMPLE_RATE) / FFT_SIZE
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

function windowFunction(index) {
  if (selectedWindow === WINDOW_FUNCTIONS.HANNING) {
    // Hanning window
    const alpha = 0.5
    return 0.5 * (1 - Math.cos((2 * Math.PI * index) / (FFT_SIZE - 1)))
  } else if (selectedWindow === WINDOW_FUNCTIONS.GAUSSIAN) {
    // Gaussian window
    const sigma = 0.4 // Adjust as needed
    const center = (FFT_SIZE - 1) / 2
    return Math.exp(-0.5 * ((index - center) / (sigma * center)) ** 2)
  } else {
    // Default to rectangular window
    return 1
  }
}

// // Event listener for window function selection
// document.getElementById('windowSelect').addEventListener('change', (event) => {
//   selectedWindow = event.target.value
// })

export { addToBuffer, calculateChromagram }
