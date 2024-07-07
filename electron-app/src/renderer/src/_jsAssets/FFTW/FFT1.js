// fftwWrapper.js

// Import the FFTWModule (assuming it's already loaded)
import * as FFTWModule from './FFTW.js' // Adjust the path accordingly

// Create the FFTW class
class FFTW {
  constructor(size) {
    this.size = size
    this.rptr = FFTWModule._malloc(size * 4 + (size + 2) * 4)
    this.cptr = this.rptr + size * 4
    this.r = new Float32Array(FFTWModule.HEAPU8.buffer, this.rptr, size)
    this.c = new Float32Array(FFTWModule.HEAPU8.buffer, this.cptr, size + 2)

    this.fplan = FFTWModule.cwrap('fftwf_plan_dft_r2c_1d', 'number', [
      'number',
      'number',
      'number',
      'number'
    ])(size, this.rptr, this.cptr, FFTWModule.FFTW_ESTIMATE)

    this.iplan = FFTWModule.cwrap('fftwf_plan_dft_c2r_1d', 'number', [
      'number',
      'number',
      'number',
      'number'
    ])(size, this.cptr, this.rptr, FFTWModule.FFTW_ESTIMATE)
  }

  forward(real) {
    this.r.set(real)
    FFTWModule.cwrap('fftwf_execute', 'void', ['number'])(this.fplan)
    return new Float32Array(FFTWModule.HEAPU8.buffer, this.cptr, this.size + 2)
  }

  inverse(cpx) {
    this.c.set(cpx)
    FFTWModule.cwrap('fftwf_execute', 'void', ['number'])(this.iplan)
    return new Float32Array(FFTWModule.HEAPU8.buffer, this.rptr, this.size)
  }

  dispose() {
    FFTWModule.cwrap('fftwf_destroy_plan', 'void', ['number'])(this.fplan)
    FFTWModule.cwrap('fftwf_destroy_plan', 'void', ['number'])(this.iplan)
    FFTWModule._free(this.rptr)
  }
}

// Create the RFFTW class (similar to FFTW)
class RFFTW {
  // ... (same as your existing RFFTW class implementation) ...
}

// Export the classes
export { FFTW, RFFTW }
