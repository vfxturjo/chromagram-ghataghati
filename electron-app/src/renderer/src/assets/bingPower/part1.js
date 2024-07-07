import FFT from 'fftw-js/src/main'

// Constants for chromagram calculation
const NUM_BINS = 12 // Number of pitch classes (C, C#, D, etc.)
const SAMPLE_RATE = 44100 // Adjust based on your audio source
const FFT_SIZE = 2048 // Size of the FFT window
const BUFFER_LEN = FFT_SIZE / 2 // Length of the circular buffer

// Circular buffer to store recent audio samples
const circularBuffer = new Float32Array(BUFFER_LEN)
let bufferIndex = 0

// Function to add samples to the circular buffer
function addToBuffer(samples) {
  samples.forEach((sample) => {
    circularBuffer[bufferIndex] = sample
    bufferIndex = (bufferIndex + 1) % BUFFER_LEN
  })
}

// Function to calculate the chromagram from the circular buffer
function calculateChromagram() {
  // Create an array to store the current window of audio samples
  const windowedSamples = new Float32Array(FFT_SIZE)

  // Copy samples from the circular buffer to the windowed array
  for (let i = 0; i < FFT_SIZE; i++) {
    const circularIndex = (bufferIndex + i) % BUFFER_LEN
    windowedSamples[i] = circularBuffer[circularIndex]
  }

  // Perform FFT on the windowed samples
  // Note: You'll need to implement FFT or use a library
  const frequencyData = performFFT(windowedSamples)

  // Map frequency data to pitch classes and accumulate energy
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

// Placeholder for FFT implementation
function performFFT(samples) {
  // Implement FFT or use a library to convert time-domain samples to frequency-domain
  return new Float32Array(FFT_SIZE) // Replace with actual frequency data
}

export { addToBuffer, calculateChromagram }
