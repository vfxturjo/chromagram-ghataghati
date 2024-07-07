<!--
Created: Wed May 29 2024 14:19:23 GMT+0600 (Bangladesh Standard Time)
Modified: Wed May 29 2024 14:44:01 GMT+0600 (Bangladesh Standard Time)
-->

# Analog-style spectrum analyzer and sliding DFT visualization using AudioWorklet 

[see here](https://codepen.io/TF3RDL/pen/MWLzPoO)

and

[see here too](https://codepen.io/TF3RDL/pen/poQJwRW)

```html
<div id="container"><canvas id="canvas"></canvas></div>
<audio id="audio" controls crossorigin></audio>
<input id="audioFileInput" type="file" accept="audio/*">
<script id="AudioProvider" type="worklet">
    class AudioProvider extends AudioWorkletProcessor {
  constructor() {
    super();
    this.dataArrays = [];
    this.bufferSize = 32768; // can handle more than 32768 samples of PCM data unlike in AnalyserNode.getFloatTimeDomainData, which is capped at 32768 samples
    this.bufferIdx = 0;
    this.currentTimeInSamples = 0;
    this.port.onmessage = (e) => {
      const audioChunks = [],
            retrievalWindowSize = e.data ? Math.min(this.bufferSize, currentFrame - this.currentTimeInSamples) : this.bufferSize,
            timeOffset = this.bufferSize-retrievalWindowSize;
      for (let channelIdx = 0; channelIdx < this.dataArrays.length; channelIdx++) {
        audioChunks[channelIdx] = [];
        for (let i = 0; i < this.dataArrays[channelIdx].length-timeOffset; i++) {
          const data = this.dataArrays[channelIdx][((this.bufferIdx+i+timeOffset) % this.bufferSize + this.bufferSize) % this.bufferSize];
          audioChunks[channelIdx][i] = data !== undefined ? data : 0;
        }
      }
      this.port.postMessage({currentChunk: audioChunks});
      this.currentTimeInSamples = currentFrame;
    };
  }
  
  process(inputs, _, _2) {
    if (inputs[0].length <= 0)
      return true;
    this.dataArrays.length = inputs[0].length;
    for (let i = 0; i < this.dataArrays.length; i++) {
      if (this.dataArrays[i] === undefined)
        this.dataArrays[i] = new Array(this.bufferSize);
      else {
        this.dataArrays[i].length = this.bufferSize;
      }
    }
    
    for (let i = 0; i < inputs[0][0].length; i++) {
      this.bufferIdx = Math.min(this.bufferIdx, this.bufferSize-1);
      for (let channelIdx = 0; channelIdx < inputs[0].length; channelIdx++) {
        this.dataArrays[channelIdx][this.bufferIdx] = inputs[0][channelIdx][i];
      }
      this.bufferIdx = ((this.bufferIdx + 1) % this.bufferSize + this.bufferSize) % this.bufferSize;
    }
    return true;
  }
}

registerProcessor('audio-provider', AudioProvider);
</script>
<script>
    class AnalogStyleAnalyzer {
        constructor(...args) {
            // initialize the sDFT coefficients
            this.calcCoeffs(args);
            this.spectrumData = [];
        }

        calcCoeffs(freqBands, order = 4, timeRes = Infinity, bandwidth = 1, sampleRate = 44100, compensateBW = true, prewarpQ = false) {
            this._coeffs = freqBands.map(x => {
                // biquad bandpass filter (cascaded biquad bandpass is not Butterworth nor Bessel, rather it is something called "critically-damped" since each filter stage shares the same every biquad coefficients)
                const K = Math.tan(Math.PI * x.ctr / sampleRate),
                    bw = Math.abs(x.hi - x.lo) * bandwidth + 1 / (timeRes / 1000),
                    qCompensationFactor = prewarpQ ? (Math.PI * x.ctr / sampleRate) / K : 1,
                    Q = x.ctr / bw * qCompensationFactor / (compensateBW ? Math.sqrt(order) : 1),
                    norm = 1 / (1 + K / Q + K * K),
                    a0 = K / Q * norm,
                    a1 = 0,
                    a2 = -a0,
                    b1 = 2 * (K * K - 1) * norm,
                    b2 = (1 - K / Q + K * K) * norm,
                    zs = [];
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
                };
            });
        }

        analyze(samples) {
            const newSpectrumData = new Array(this._coeffs.length).fill(0);
            for (const x of samples) {
                for (let i = 0; i < this._coeffs.length; i++) {
                    for (let j = 0; j < this._coeffs[i].zs.length; j++) {
                        const input = j <= 0 ? x : this._coeffs[i].zs[j - 1].out;
                        this._coeffs[i].zs[j].out = input * this._coeffs[i].a0 + this._coeffs[i].zs[j].z1;
                        this._coeffs[i].zs[j].z1 = input * this._coeffs[i].a1 + this._coeffs[i].zs[j].z2 - this._coeffs[i].b1 * this._coeffs[i].zs[j].out;
                        this._coeffs[i].zs[j].z2 = input * this._coeffs[i].a2 - this._coeffs[i].b2 * this._coeffs[i].zs[j].out;
                    }
                    newSpectrumData[i] = Math.max(newSpectrumData[i], Math.abs(this._coeffs[i].zs[this._coeffs[i].zs.length - 1].out));
                }
            }
            this.spectrumData = newSpectrumData.map(x => x / 2);
        }
    }
</script>
<script>
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
            this.calcCoeffs(args);
            this.spectrumData = [];
        }

        calcCoeffs(freqBands, order = 4, timeRes = 600, bandwidth = 1, sampleRate = 44100, compensateBW = true) {
            // calcCoeffs() can be called anywhere else to re-initialize sliding DFT after changes in frequency band distributions and note that x and y are used instead of real and imaginary since vector rotation is the equivalent of the complex one
            this._coeffs = [];
            freqBands.map((x, i) => {
                // rX and rY are calculated in advance here since calculating sin and cos functions are pretty slow af
                this._coeffs[i] = {
                    rX: Math.cos(x.ctr * Math.PI / sampleRate * 2),
                    rY: Math.sin(x.ctr * Math.PI / sampleRate * 2),
                    decay: Math.E ** ((-Math.abs(x.hi - x.lo) * Math.PI * bandwidth / sampleRate - 1 / (timeRes * sampleRate / (Math.PI * 1000))) * (compensateBW ? Math.sqrt(order) : 1)),
                    coeffs: []
                };
                for (let j = 0; j < order; j++) {
                    this._coeffs[i].coeffs[j] = {
                        x: 0,
                        y: 0
                    };
                }
            });
        }

        analyze(dataArray) {
            const newSpectrumData = new Array(this._coeffs.length).fill(0);
            for (const x of dataArray) {
                for (let i = 0; i < this._coeffs.length; i++) {
                    for (let j = 0; j < this._coeffs[i].coeffs.length; j++) {
                        const input = j <= 0 ? {
                                x: x,
                                y: 0,
                            } : this._coeffs[i].coeffs[j - 1],
                            outX = (this._coeffs[i].coeffs[j].x * this._coeffs[i].rX - this._coeffs[i].coeffs[j].y * this._coeffs[i].rY) * this._coeffs[i].decay + input.x * (1 - this._coeffs[i].decay),
                            outY = (this._coeffs[i].coeffs[j].x * this._coeffs[i].rY + this._coeffs[i].coeffs[j].y * this._coeffs[i].rX) * this._coeffs[i].decay + input.y * (1 - this._coeffs[i].decay);

                        this._coeffs[i].coeffs[j].x = outX;
                        this._coeffs[i].coeffs[j].y = outY;
                    }
                    newSpectrumData[i] = Math.max(newSpectrumData[i],
                        this._coeffs[i].coeffs[this._coeffs[i].coeffs.length - 1].x ** 2 +
                        this._coeffs[i].coeffs[this._coeffs[i].coeffs.length - 1].y ** 2);
                }
            }
            this.spectrumData = newSpectrumData.map((x) => Math.sqrt(x));
        }
    }
</script>
<script>
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
    class VQsDFT {
        constructor(...args) {
            this.calcCoeffs(args);
            this.spectrumData = [];
        }

        calcCoeffs(freqBands, window = [1, 0.5], timeRes = 600, bandwidth = 1, bufferSize = 44100, sampleRate = 44100, useNC = false) {
            this._coeffs = freqBands.map(x => {
                const fiddles = [],
                    twiddles = [],
                    resonCoeffs = [],
                    coeffs1 = [],
                    coeffs2 = [],
                    coeffs3 = [],
                    coeffs4 = [],
                    coeffs5 = [],
                    gains = [],
                    period = Math.trunc(Math.min(bufferSize, sampleRate / (bandwidth * Math.abs(x.hi - x.lo) + 1 / (timeRes / 1000)))), // N must be an integer, but K doesn't have to be
                    minIdx = useNC ? 0 : -window.length + 1,
                    maxIdx = useNC ? 2 : window.length;
                // this below is needed since we have to apply a frequency-domain window function
                for (let i = minIdx; i < maxIdx; i++) {
                    const amplitude = useNC ? 1 : window[Math.abs(i)] * (-(Math.abs(i) % 2) * 2 + 1),
                        k = x.ctr * period / sampleRate + i - useNC / 2,
                        fid = -2 * Math.PI * k,
                        twid = 2 * Math.PI * k / period,
                        reson = 2 * Math.cos(2 * Math.PI * k / period);
                    fiddles.push({
                        x: Math.cos(fid),
                        y: Math.sin(fid)
                    });
                    twiddles.push({
                        x: Math.cos(twid),
                        y: Math.sin(twid)
                    });
                    resonCoeffs.push(reson);
                    coeffs1.push({
                        x: 0,
                        y: 0
                    });
                    coeffs2.push({
                        x: 0,
                        y: 0
                    });
                    coeffs3.push({
                        x: 0,
                        y: 0
                    });
                    coeffs4.push({
                        x: 0,
                        y: 0
                    });
                    coeffs5.push({
                        x: 0,
                        y: 0
                    });
                    gains.push(amplitude);
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
                };
            });
            this._buffer = new Array(bufferSize + 1).fill(0);
            this._bufferIdx = this._buffer.length - 1; // this is required for circular buffer
        }

        analyze(samples) {
            this.spectrumData = new Array(this._coeffs.length).fill(0);
            for (const sample of samples) {
                // Admittedly slow linear buffer
                /*
                this._buffer.push(sample);
                this._buffer.shift();
                */
                // Circular buffer
                this._bufferIdx = ((this._bufferIdx + 1) % this._buffer.length + this._buffer.length) % this._buffer.length;
                this._buffer[this._bufferIdx] = sample;
                for (let i = 0; i < this._coeffs.length; i++) {
                    const coeff = this._coeffs[i],
                        kernelLength = coeff.coeffs1.length,
                        /*oldest = this._buffer.length-coeff.period-1,
                        latest = this._buffer.length-1,*/
                        oldest = ((this._bufferIdx - coeff.period) % this._buffer.length + this._buffer.length) % this._buffer.length,
                        latest = this._bufferIdx,
                        sum = {
                            x: 0,
                            y: 0
                        };
                    for (let j = 0; j < kernelLength; j++) {
                        const fiddle = coeff.fiddles[j],
                            twiddle = coeff.twiddles[j],
                            // Comb stage
                            combX = this._buffer[latest] * fiddle.x - this._buffer[oldest],
                            combY = this._buffer[latest] * fiddle.y

                        // Second stage
                        coeff.coeffs1[j].x = combX * twiddle.x - combY * twiddle.y - coeff.coeffs2[j].x;
                        coeff.coeffs1[j].y = combX * twiddle.y + combY * twiddle.x - coeff.coeffs2[j].y;

                        coeff.coeffs2[j].x = combX;
                        coeff.coeffs2[j].y = combY;

                        // Real resonator
                        coeff.coeffs3[j].x = coeff.coeffs1[j].x + coeff.resonCoeffs[j] * coeff.coeffs4[j].x - coeff.coeffs5[j].x;
                        coeff.coeffs3[j].y = coeff.coeffs1[j].y + coeff.resonCoeffs[j] * coeff.coeffs4[j].y - coeff.coeffs5[j].y;

                        coeff.coeffs5[j].x = coeff.coeffs4[j].x;
                        coeff.coeffs5[j].y = coeff.coeffs4[j].y;

                        coeff.coeffs4[j].x = coeff.coeffs3[j].x;
                        coeff.coeffs4[j].y = coeff.coeffs3[j].y;

                        sum.x += coeff.coeffs3[j].x * coeff.gains[j] / coeff.period;
                        sum.y += coeff.coeffs3[j].y * coeff.gains[j] / coeff.period;
                    }
                    const period = coeff.period
                    this.spectrumData[i] = Math.max(this.spectrumData[i], coeff.nc ? -(coeff.coeffs3[0].x / period * coeff.coeffs3[1].x / period) - (coeff.coeffs3[0].y / period * coeff.coeffs3[1].y / period) : sum.x ** 2 + sum.y ** 2);
                }
            }
            this.spectrumData = this.spectrumData.map(x => Math.sqrt(x));
        }
    }
</script>
<script>
    function map(x, min, max, targetMin, targetMax) {
        return (x - min) / (max - min) * (targetMax - targetMin) + targetMin;
    }

    function clamp(x, min, max) {
        return Math.min(Math.max(x, min), max);
    }

    function idxWrapOver(x, length) {
        return (x % length + length) % length;
    }
    // Hz and FFT bin conversion
    function hertzToFFTBin(x, y = 'round', bufferSize = 4096, sampleRate = 44100) {
        const bin = x * bufferSize / sampleRate;
        let func = y;

        if (!['floor', 'ceil', 'trunc'].includes(func))
            func = 'round'; // always use round if you specify an invalid/undefined value

        return Math[func](bin);
    }

    function fftBinToHertz(x, bufferSize = 4096, sampleRate = 44100) {
        return x * sampleRate / bufferSize;
    }

    // Calculate the FFT
    function calcFFT(input) {
        let fft = input.map(x => x);
        let fft2 = input.map(x => x);
        transform(fft, fft2);
        let output = new Array(Math.round(fft.length / 2)).fill(0);
        for (let i = 0; i < output.length; i++) {
            output[i] = Math.hypot(fft[i], fft2[i]) / (fft.length);
        }
        return output;
    }

    function calcComplexFFT(input) {
        let fft = input.map(x => x);
        let fft2 = input.map(x => x);
        transform(fft, fft2);
        return input.map((_, i, arr) => {
            return {
                re: fft[i] / (arr.length / 2),
                im: fft2[i] / (arr.length / 2),
                magnitude: Math.hypot(fft[i], fft2[i]) / (arr.length / 2),
                phase: Math.atan2(fft2[i], fft[i])
            };
        });
    }

    function calcComplexInputFFT(real, imag) {
        if (real.length !== imag.length)
            return [];
        const fft1 = real.map(x => x),
            fft2 = imag.map(x => x);
        transform(fft1, fft2);
        return real.map((_, i, arr) => {
            return {
                re: fft1[i] / arr.length,
                im: fft2[i] / arr.length,
                magnitude: Math.hypot(fft1[i], fft2[i]) / arr.length,
                phase: Math.atan2(fft2[i], fft1[i])
            }
        });
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
        const n = real.length;
        if (n != imag.length)
            throw "Mismatched lengths";
        if (n <= 0)
            return;
        else if ((2 ** Math.trunc(Math.log2(n))) === n) // Is power of 2
            transformRadix2(real, imag);
        else // More complicated algorithm for arbitrary sizes
            transformBluestein(real, imag);
    }

    /* 
     * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
     * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
     */
    function inverseTransform(real, imag) {
        transform(imag, real);
    }

    /* 
     * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
     * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
     */
    function transformRadix2(real, imag) {
        // Length variables
        const n = real.length;
        if (n != imag.length)
            throw "Mismatched lengths";
        if (n <= 1) // Trivial transform
            return;
        const logN = Math.log2(n);
        if ((2 ** Math.trunc(logN)) !== n)
            throw "Length is not a power of 2";

        // Trigonometric tables
        let cosTable = new Array(n / 2);
        let sinTable = new Array(n / 2);
        for (let i = 0; i < n / 2; i++) {
            cosTable[i] = Math.cos(2 * Math.PI * i / n);
            sinTable[i] = Math.sin(2 * Math.PI * i / n);
        }

        // Bit-reversed addressing permutation
        for (let i = 0; i < n; i++) {
            let j = reverseBits(i, logN);
            if (j > i) {
                let temp = real[i];
                real[i] = real[j];
                real[j] = temp;
                temp = imag[i];
                imag[i] = imag[j];
                imag[j] = temp;
            }
        }

        // Cooley-Tukey decimation-in-time radix-2 FFT
        for (let size = 2; size <= n; size *= 2) {
            let halfsize = size / 2;
            let tablestep = n / size;
            for (let i = 0; i < n; i += size) {
                for (let j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
                    const l = j + halfsize;
                    const tpre = real[l] * cosTable[k] + imag[l] * sinTable[k];
                    const tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
                    real[l] = real[j] - tpre;
                    imag[l] = imag[j] - tpim;
                    real[j] += tpre;
                    imag[j] += tpim;
                }
            }
        }

        // Returns the integer whose value is the reverse of the lowest 'bits' bits of the integer 'x'.
        function reverseBits(x, bits) {
            let y = 0;
            for (let i = 0; i < bits; i++) {
                y = (y << 1) | (x & 1);
                x >>>= 1;
            }
            return y;
        }
    }

    /* 
     * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
     * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
     * Uses Bluestein's chirp z-transform algorithm.
     */
    function transformBluestein(real, imag) {
        // Find a power-of-2 convolution length m such that m >= n * 2 + 1
        const n = real.length;
        if (n != imag.length)
            throw "Mismatched lengths";
        const m = 2 ** Math.trunc(Math.log2(n * 2) + 1);

        // Trignometric tables
        let cosTable = new Array(n);
        let sinTable = new Array(n);
        for (let i = 0; i < n; i++) {
            let j = i * i % (n * 2); // This is more accurate than j = i * i
            cosTable[i] = Math.cos(Math.PI * j / n);
            sinTable[i] = Math.sin(Math.PI * j / n);
        }

        // Temporary vectors and preprocessing
        let areal = newArrayOfZeros(m);
        let aimag = newArrayOfZeros(m);
        for (let i = 0; i < n; i++) {
            areal[i] = real[i] * cosTable[i] + imag[i] * sinTable[i];
            aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];
        }
        let breal = newArrayOfZeros(m);
        let bimag = newArrayOfZeros(m);
        breal[0] = cosTable[0];
        bimag[0] = sinTable[0];
        for (let i = 1; i < n; i++) {
            breal[i] = breal[m - i] = cosTable[i];
            bimag[i] = bimag[m - i] = sinTable[i];
        }

        // Convolution
        let creal = new Array(m);
        let cimag = new Array(m);
        convolveComplex(areal, aimag, breal, bimag, creal, cimag);

        // Postprocessing
        for (let i = 0; i < n; i++) {
            real[i] = creal[i] * cosTable[i] + cimag[i] * sinTable[i];
            imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];
        }
    }

    /* 
     * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
     */
    function convolveReal(x, y, out) {
        const n = x.length;
        if (n != y.length || n != out.length)
            throw "Mismatched lengths";
        convolveComplex(x, newArrayOfZeros(n), y, newArrayOfZeros(n), out, newArrayOfZeros(n));
    }

    /* 
     * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
     */
    function convolveComplex(xreal, ximag, yreal, yimag, outreal, outimag) {
        const n = xreal.length;
        if (n != ximag.length || n != yreal.length || n != yimag.length ||
            n != outreal.length || n != outimag.length)
            throw "Mismatched lengths";

        xreal = xreal.slice();
        ximag = ximag.slice();
        yreal = yreal.slice();
        yimag = yimag.slice();
        transform(xreal, ximag);
        transform(yreal, yimag);

        for (let i = 0; i < n; i++) {
            const temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
            ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
            xreal[i] = temp;
        }
        inverseTransform(xreal, ximag);

        for (let i = 0; i < n; i++) { // Scaling (because this FFT implementation omits it)
            outreal[i] = xreal[i] / n;
            outimag[i] = ximag[i] / n;
        }
    }

    function newArrayOfZeros(n) {
        let result = new Array(n).fill(0);
        return result;
    }
</script>
```

```css
body {
    margin: 0;
    overflow: hidden;
}

audio {
    display: inline-block;
    width: 100%;
    height: 40px;
}

canvas {
    display: block;
    width: 100%;
}

#container {
    height: calc(100vh - 40px);
}

#upload {
    display: none;
}
```

```javascript
// necessary parts for audio context and audio elements respectively
const audioCtx = new AudioContext();
const audioPlayer = document.getElementById('audio');
const localAudioElement = document.getElementById('audioFileInput');
localAudioElement.addEventListener('change', loadLocalFile);
// canvas is for displaying visuals
const canvas = document.getElementById('canvas'),
    ctx = canvas.getContext('2d'),
    container = document.getElementById('container');
const audioSource = audioCtx.createMediaElementSource(audioPlayer);
const analyser = audioCtx.createAnalyser();
analyser.fftSize = 32768; // maxes out FFT size
const dataArray = new Float32Array(analyser.fftSize);
// variables
const currentSpectrum = [],
    peaks = [],
    peakHolds = [],
    peakAccels = [];
const delay = audioCtx.createDelay();
audioSource.connect(delay);
delay.connect(audioCtx.destination);
//audioSource.connect(audioCtx.destination);
audioSource.connect(analyser);
let audioProvider,
    currentSampleRate = audioCtx.sampleRate,
    freqBands = [];
const analogStyleAnalyser = new AnalogStyleAnalyzer([]),
    swift = new SWIFT([]),
    sdft = new VQsDFT([]);
const customDSPSource = document.getElementById('AudioProvider'),
    dspSourceBlob = new Blob([customDSPSource.innerText], {
        type: 'application/javascript'
    }),
    dspSourceUrl = URL.createObjectURL(dspSourceBlob);

const visualizerSettings = {
        //fftSize: 1152,
        freqDist: 'octaves',
        numBands: 50, // similar to WMP's Bars visualization when number of bands are at maximum possible
        minFreq: 20,
        maxFreq: 20000,
        fscale: 'logarithmic',
        hzLinearFactor: 0,
        minNote: 4,
        maxNote: 124,
        noteTuning: 1000, // setting it to 1kHz does automatically makes octave bands compliant with ANSI S1.11-2004 standard when comes to one-third octave band center frequencies right?
        octaves: 6, // defaults to something similar to Spectroscope visualization in WaveLab
        detune: 0,
        analysisAlgorithm: 'analog',
        bandwidth: 1,
        order: 1,
        prewarpQ: true,
        compensateBW: true,
        windowFunction: '1, 0.5',
        customWindow: '1',
        useNC: false,
        timeRes: 100,
        maxTimeRes: 1000,
        constantQ: true,
        resetCoeffs: recalcCoeffs,
        resetAverages: resetSmoothedValues,
        useAccurateSmoothing: true,
        smoothingTimeConstant: 90, // default value is approximately the main bar of audio visualizer thing in Geometry Dash 2.2
        useAverageSmoothing: false,
        peakDecay: 0,
        peakHold: 60,
        useActualPeak: false,
        fadingPeaks: true, // this effect is used on peak hold part of Audio Visualizer effect on GD 2.2
        // minDecibels and maxDecibels defaults to -60...+6 to match foobar2000's built-in Spectrum visualization
        minDecibels: -60,
        maxDecibels: 6,
        useDecibels: true,
        gamma: 1,
        useAbsolute: true,
        showLabels: true,
        showLabelsY: true,
        amplitudeLabelInterval: 10,
        labelTuning: 440,
        showDC: true,
        showNyquist: true,
        mirrorLabels: true,
        diffLabels: false,
        labelMode: 'decade',
        freeze: false,
        useGradient: true,
        darkMode: false,
        showPeaks: true,
        resetBoth: resetBoth,
        useIncorrectWay: false, // when enabled, it uses the wrong way of getting samples
        fftSize: 576 // default is the buffer length of PCM data on Winamp's visualization system
        //compensateDelay: true
    },
    placeholderListData = {
        'Option 1': 'one',
        'Option 2': 'two'
    },
    loader = {
        url: '',
        load: function() {
            audioPlayer.src = this.url;
            audioPlayer.play();
        },
        loadLocal: function() {
            localAudioElement.click();
        },
        toggleFullscreen: _ => {
            if (document.fullscreenElement === canvas)
                document.exitFullscreen();
            else
                canvas.requestFullscreen();
        }
    };
// dat.GUI for quick customization
let gui = new dat.GUI();
gui.add(loader, 'url').name('URL');
gui.add(loader, 'load').name('Load');
gui.add(loader, 'loadLocal').name('Load from local device');
let settings = gui.addFolder('Visualization settings');
// FFT size can be non-power of 2 because we use the FFT library that supports non-power of two data length
//settings.add(visualizerSettings, 'fftSize', 32, 32768, 1).name('FFT size');
// The additional parameters goes here
// another parameters at the end
const freqDistFolder = settings.addFolder('Frequency distribution');
freqDistFolder.add(visualizerSettings, 'freqDist', {
    'Frequency bands': 'freqs',
    'Octave bands': 'octaves'
}).name('Frequency band distribution').onChange(recalcCoeffs);
// up to 192kHz sample rate is supported for full-range visualization
freqDistFolder.add(visualizerSettings, 'minFreq', 0, 96000).name('Minimum frequency').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'maxFreq', 0, 96000).name('Maximum frequency').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'minNote', 0, 128).name('Minimum note').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'maxNote', 0, 128).name('Maximum note').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'noteTuning', 0, 96000).name('Octave bands tuning (nearest note = tuning frequency in Hz)').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'detune', -24, 24).name('Detune').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'octaves', 1, 48).name('Bands per octave').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'numBands', 2, 512, 1).name('Number of bands').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'fscale', {
    'Bark': 'bark',
    'ERB': 'erb',
    'Cams': 'cam',
    'Mel (AIMP)': 'mel',
    'Linear': 'linear',
    'Logarithmic': 'logarithmic',
    'Hyperbolic sine': 'sinh',
    'Shifted logarithmic': 'shifted log',
    'Nth root': 'nth root',
    'Negative exponential': 'negative exponential',
    'Adjustable Bark': 'adjustable bark',
    'Period': 'period'
}).name('Frequency scale').onChange(recalcCoeffs);
freqDistFolder.add(visualizerSettings, 'hzLinearFactor', 0, 100).name('Hz linear factor').onChange(recalcCoeffs);
const transformFolder = settings.addFolder('Transform algorithm');
transformFolder.add(visualizerSettings, 'analysisAlgorithm', {
    'Analog-style analyzer': 'analog',
    'Sliding windowed infinite Fourier transform': 'swift',
    'Variable-Q sliding DFT': 'sdft'
}).name('Analysis algorithm').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'bandwidth', 0, 64).name('Bandwidth').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'order', 1, 8, 1).name('Filter order').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'windowFunction', {
    'Rectangular': '1',
    'Hann': '1, 0.5',
    'Hamming': '1, 0.4259434938430786',
    'Blackman': '1, 0.595257580280304, 0.0952545627951622',
    'Nuttall': '1, 0.6850073933601379, 0.20272639393806458, 0.017719272524118423',
    'Flat top': '1, 0.966312825679779, 0.6430955529212952, 0.19387830793857574, 0.016120079904794693',
    'Custom': 'custom'
}).name('Window function').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'customWindow').name('Custom frequency-domain windowing coefficients').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'useNC').name('Use NC method (VQ-sDFT only)').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'prewarpQ').name('Use prewarped Q (analog-style analyzer only)').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'compensateBW').name('Compensate bandwidth for narrowing on higher order filters (IIR filter banks only)').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'timeRes', 0, 2000).name('Time resoluion').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'constantQ').name('Use constant-Q instead of variable-Q').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'maxTimeRes', 0, 8000).name('Maximum time resoluion').onChange(recalcCoeffs);
transformFolder.add(visualizerSettings, 'resetCoeffs').name('Reset coefficients');
const peakFolder = settings.addFolder('Time averaging and peak decay settings');
peakFolder.add(visualizerSettings, 'useAccurateSmoothing').name('Apply time smoothing during processing');
peakFolder.add(visualizerSettings, 'smoothingTimeConstant', 0, 100).name('Smoothing time constant');
peakFolder.add(visualizerSettings, 'useAverageSmoothing').name('Use exponential average instead of peak decay');
peakFolder.add(visualizerSettings, 'peakHold', 0, 120).name('Peak hold time');
peakFolder.add(visualizerSettings, 'peakDecay', 0, 100).name('Peak fall rate');
peakFolder.add(visualizerSettings, 'useActualPeak').name('Use actual peak');
peakFolder.add(visualizerSettings, 'resetAverages').name('Reset smoothed values and peaks');
const amplitudeFolder = settings.addFolder('Amplitude');
amplitudeFolder.add(visualizerSettings, 'useDecibels').name('Use logarithmic amplitude/decibel scale');
amplitudeFolder.add(visualizerSettings, 'useAbsolute').name('Use absolute value');
amplitudeFolder.add(visualizerSettings, 'gamma', 0.5, 10).name('Gamma');
amplitudeFolder.add(visualizerSettings, 'minDecibels', -120, 6).name('Lower amplitude range');
amplitudeFolder.add(visualizerSettings, 'maxDecibels', -120, 6).name('Higher amplitude range');
const labelFolder = settings.addFolder('Labels and grids');
labelFolder.add(visualizerSettings, 'showLabels').name('Show horizontal-axis labels');
labelFolder.add(visualizerSettings, 'showLabelsY').name('Show vertical-axis labels');
labelFolder.add(visualizerSettings, 'amplitudeLabelInterval', 0.5, 48).name('dB label interval');
labelFolder.add(visualizerSettings, 'showDC').name('Show DC label');
labelFolder.add(visualizerSettings, 'showNyquist').name('Show Nyquist frequency label');
labelFolder.add(visualizerSettings, 'mirrorLabels').name('Mirror Y-axis labels');
labelFolder.add(visualizerSettings, 'diffLabels').name('Use difference coloring for labels');
labelFolder.add(visualizerSettings, 'labelMode', {
    'Decades': 'decade',
    'Octaves': 'octave',
    'Notes': 'note',
    'Automatic': 'auto'
}).name('Frequency label mode');
labelFolder.add(visualizerSettings, 'labelTuning', 0, 96000).name('Note labels tuning (nearest note = tuning frequency in Hz)');
settings.add(visualizerSettings, 'showPeaks').name('Show peaks');
settings.add(visualizerSettings, 'fadingPeaks').name('Enable peak fading effect');
settings.add(visualizerSettings, 'freeze').name('Freeze analyzer');
settings.add(visualizerSettings, 'useGradient').name('Use color gradient');
settings.add(visualizerSettings, 'darkMode').name('Dark mode');
settings.add(visualizerSettings, 'useIncorrectWay').name('Use getFloatTimeDomainData instead of AudioWorklet');
settings.add(visualizerSettings, 'fftSize', 32, 32768, 1).name('getFloatTimeDomainData buffer length (samples)');
settings.add(visualizerSettings, 'resetBoth').name('Reset both coefficients and smoothing');
//settings.add(visualizerSettings, 'compensateDelay').name('Compensate for delay');
gui.add(loader, 'toggleFullscreen').name('Toggle fullscreen mode');

function resetBoth() {
    resetSmoothedValues();
    recalcCoeffs();
}

function resetSmoothedValues() {
    updateSpectrumVisualization([]);
}

function recalcCoeffs() {
    switch (visualizerSettings.freqDist) {
        case 'octaves':
            freqBands = generateOctaveBands(visualizerSettings.octaves, visualizerSettings.minNote, visualizerSettings.maxNote, visualizerSettings.detune, visualizerSettings.noteTuning);
            break;
        default:
            freqBands = generateFreqBands(visualizerSettings.numBands, visualizerSettings.minFreq, visualizerSettings.maxFreq, visualizerSettings.fscale, visualizerSettings.hzLinearFactor / 100);
    }

    const windowingKernel = parseList(visualizerSettings.windowFunction === 'custom' ? visualizerSettings.customWindow : visualizerSettings.windowFunction),
        timeRes = visualizerSettings.constantQ ? Infinity : visualizerSettings.timeRes,
        iirArgs = [freqBands, visualizerSettings.order, timeRes, visualizerSettings.bandwidth, audioCtx.sampleRate, visualizerSettings.compensateBW, visualizerSettings.prewarpQ],
        firArgs = [freqBands, windowingKernel, timeRes, visualizerSettings.bandwidth, Math.round(audioCtx.sampleRate * visualizerSettings.maxTimeRes / 1000), audioCtx.sampleRate, visualizerSettings.useNC];
    analogStyleAnalyser.calcCoeffs([]);
    swift.calcCoeffs([]);
    sdft.calcCoeffs([]);
    switch (visualizerSettings.analysisAlgorithm) {
        case 'analog':
            analogStyleAnalyser.calcCoeffs(...iirArgs);
        case 'swift':
            swift.calcCoeffs(...iirArgs);
        default:
            sdft.calcCoeffs(...firArgs);
    }
}
recalcCoeffs();

function resizeCanvas() {
    const scale = devicePixelRatio,
        isFullscreen = document.fullscreenElement === canvas;
    canvas.width = (isFullscreen ? innerWidth : container.clientWidth) * scale;
    canvas.height = (isFullscreen ? innerHeight : container.clientHeight) * scale;
}

addEventListener('click', () => {
    if (audioCtx.state == 'suspended')
        audioCtx.resume();
});
addEventListener('resize', resizeCanvas);
resizeCanvas();

function loadLocalFile(event) {
    const file = event.target.files[0],
        reader = new FileReader();
    reader.onload = (e) => {
        audioPlayer.src = e.target.result;
        audioPlayer.play();
    };

    reader.readAsDataURL(file);
}

//visualize();
audioCtx.audioWorklet.addModule(dspSourceUrl).then(() => {
    //let messageCounter = 0;
    audioProvider = new AudioWorkletNode(audioCtx, 'audio-provider');
    audioSource.connect(audioProvider);
    audioProvider.port.postMessage(0);
    audioProvider.port.onmessage = (e) => {
        if (!visualizerSettings.freeze && !visualizerSettings.useIncorrectWay)
            analyzeChunk(e.data.currentChunk);
        audioProvider.port.postMessage(1);
        //if (messageCounter < 1) {
        //  console.log(e.data.currentChunk);
        //}
        //messageCounter++;
    };
    audioProvider.onprocessorerror = (e) => {
        console.log(e.message);
    }
    visualize();
}).catch((e) => {
    console.log(e.message);
});

function analyzeChunk(data) {
    const dataset = [];
    let retrievalLength = 0;
    for (const x of data) {
        retrievalLength = Math.max(retrievalLength, x.length);
    }
    for (let i = 0; i < retrievalLength; i++) {
        let sum = 0;
        for (let channelIdx = 0; channelIdx < data.length; channelIdx++) {
            sum += data[channelIdx][i];
        }
        dataset[i] = sum / data.length;
        if (visualizerSettings.useAccurateSmoothing) {
            getKindofsDFT().analyze([isFinite(dataset[i]) ? dataset[i] : 0]);
            updateSpectrumVisualization(getKindofsDFT().spectrumData, true);
        }
    }
    if (dataset.length > 0 && !visualizerSettings.useAccurateSmoothing)
        getKindofsDFT().analyze(dataset.map(x => isFinite(x) ? x : 0));
}

function getKindofsDFT() {
    switch (visualizerSettings.analysisAlgorithm) {
        case 'analog':
            return analogStyleAnalyser;
        case 'swift':
            return swift;
        default:
            return sdft;
    }
}

function visualize() {
    delay.delayTime.value = 0 //(visualizerSettings.fftSize / audioCtx.sampleRate) * visualizerSettings.compensateDelay;
    if (!visualizerSettings.freeze) {
        // we use getFloatTimeDomainData (which is PCM data that is gathered, just like vis_stream::get_chunk_absolute() in foobar2000 SDK)
        if (visualizerSettings.useIncorrectWay) {
            analyser.getFloatTimeDomainData(dataArray);
            const fftData = [];
            for (let i = 0; i < visualizerSettings.fftSize; i++) {
                fftData[i] = dataArray[i + analyser.fftSize - visualizerSettings.fftSize];
            }
            analyzeChunk([fftData]);
        }
        /*
        const spectrumData = getKindofsDFT().spectrumData;
        */
        if (currentSampleRate !== audioCtx.sampleRate)
            recalcCoeffs();
        /*
        currentSpectrum.length = spectrumData.length;
        for (let i = 0; i < spectrumData.length; i++) {
          currentSpectrum[i] = spectrumData[i];
        }
        */
        if (!visualizerSettings.useAccurateSmoothing)
            updateSpectrumVisualization(getKindofsDFT().spectrumData);
    }
    const fgColor = visualizerSettings.darkMode ? (visualizerSettings.useGradient ? '#c0c0c0' : '#fff') : '#000',
        bgColor = visualizerSettings.darkMode ? (visualizerSettings.useGradient ? '#202020' : '#000') : '#fff';
    let grad = fgColor;
    if (visualizerSettings.useGradient) {
        grad = ctx.createLinearGradient(0, 0, 0, canvas.height);
        // color gradient derived from foobar2000
        grad.addColorStop(0, visualizerSettings.darkMode ? '#569cd6' : 'rgb(0, 102, 204)');
        grad.addColorStop(1, visualizerSettings.darkMode ? '#c0c0c0' : '#000');
    }
    ctx.globalCompositeOperation = 'source-over';
    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = grad;
    ctx.strokeStyle = grad;
    for (let i = 0; i < currentSpectrum.length; i++) {
        ctx.fillRect(i * canvas.width / currentSpectrum.length + 1, canvas.height, canvas.width / currentSpectrum.length - 2, -map(ascale(currentSpectrum[i] * 2), 0, 1, 0, canvas.height));
    }
    ctx.fillStyle = fgColor;
    ctx.strokeStyle = fgColor;
    if (visualizerSettings.showPeaks) {
        for (let i = 0; i < peaks.length; i++) {
            ctx.globalAlpha = visualizerSettings.fadingPeaks ? peakHolds[i] / (visualizerSettings.peakHold * (visualizerSettings.useAccurateSmoothing ? audioCtx.sampleRate / 60 : 1)) : 1;
            ctx.fillRect(i * canvas.width / peaks.length + 1, map(ascale(peaks[i] * 2), 0, 1, canvas.height, 0), canvas.width / peaks.length - 2, 2);
        }
    }
    /*
    for (let i = 0; i < 24; i++) {
      let sum = 0;
      for (let j = 0; j < currentSpectrum.length/24; j++) {
        sum += currentSpectrum[i+j*24] ** 2;
      }
      ctx.fillRect(i*canvas.width/24+1, canvas.height, canvas.width/24 - 2, -Math.min(Math.sqrt(sum)*canvas.height*Math.SQRT2, canvas.height/2));
    }
    */
    ctx.globalAlpha = 1;
    ctx.globalCompositeOperation = visualizerSettings.diffLabels ? 'difference' : 'source-over';
    ctx.fillStyle = visualizerSettings.diffLabels ? '#fff' : fgColor;
    ctx.strokeStyle = visualizerSettings.diffLabels ? '#fff' : fgColor;
    // label part
    ctx.font = `${Math.trunc(10*devicePixelRatio)}px sans-serif`;
    ctx.textAlign = 'start';
    // Frequency label part
    if (visualizerSettings.showLabels || visualizerSettings.showDC || visualizerSettings.showNyquist) {
        ctx.globalAlpha = 0.5;
        ctx.setLineDash([]);

        const freqLabels = [],
            isNote = visualizerSettings.labelMode === 'note',
            notes = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'],
            minLabelRange = freqBands.length > 0 ? freqBands[0].ctr : 0,
            maxLabelRange = freqBands.length > 0 ? freqBands[freqBands.length - 1].ctr : 0,
            labelScale = visualizerSettings.freqDist === 'octaves' ? 'log' : visualizerSettings.fscale,
            hzLinearFactor = visualizerSettings.hzLinearFactor / 100;

        let freqsTable;
        switch (visualizerSettings.labelMode) {
            case 'decade':
                freqsTable = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000];
                break;
            case 'octave':
                freqsTable = [31, 63.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000];
                break;
            case 'note':
                freqsTable = generateOctaveBands(12, 0, 132, 0, visualizerSettings.labelTuning).map(x => x.ctr);
                break;
            default:
                freqsTable = freqBands.map(x => x.ctr);
        }
        if (visualizerSettings.showLabels)
            freqLabels.push(...freqsTable);
        if (visualizerSettings.showDC)
            freqLabels.push(0);
        if (visualizerSettings.showNyquist)
            freqLabels.push(audioCtx.sampleRate / 2);

        freqLabels.map(x => {
            const note = isFinite(Math.log2(x)) ? notes[idxWrapOver(Math.round(Math.log2(x) * 12), notes.length)] : 'DC',
                isSharp = note.includes('#'),
                isC = note === 'C';

            ctx.globalAlpha = isNote ? (isSharp ? 0.2 : isC ? 0.8 : 0.5) : 0.5;
            const label = x === audioCtx.sampleRate / 2 && visualizerSettings.showNyquist ? 'Nyquist' : isNote || x === 0 ? `${note}${isC ? Math.trunc(Math.log2(x)-4) : ''}` : (x >= 1000) ? `${x / 1000}kHz` : `${x}Hz`,
                posX = map(fscale(x, labelScale, hzLinearFactor), fscale(minLabelRange, labelScale, hzLinearFactor), fscale(maxLabelRange, labelScale, hzLinearFactor), canvas.width / freqBands.length / 2, canvas.width - canvas.width / freqBands.length / 2);
            ctx.beginPath();
            ctx.lineTo(posX, canvas.height);
            ctx.lineTo(posX, 0);
            ctx.stroke();
            ctx.globalAlpha = 1;
            ctx.fillText(label, posX, canvas.height);
        });
        ctx.setLineDash([]);
        ctx.globalAlpha = 1;
        ctx.textAlign = 'start';
    }
    // Amplitude/decibel label part
    if (visualizerSettings.showLabelsY) {
        const dBLabelData = [-Infinity],
            mindB = Math.min(visualizerSettings.minDecibels, visualizerSettings.maxDecibels),
            maxdB = Math.max(visualizerSettings.minDecibels, visualizerSettings.maxDecibels),
            minLabelIdx = Math.round(mindB / visualizerSettings.amplitudeLabelInterval),
            maxLabelIdx = Math.round(maxdB / visualizerSettings.amplitudeLabelInterval);

        if (isFinite(minLabelIdx) && isFinite(maxLabelIdx)) {
            for (let i = maxLabelIdx; i >= minLabelIdx; i--) {
                dBLabelData.push(i * visualizerSettings.amplitudeLabelInterval);
            }
        }

        ctx.globalAlpha = 0.5;
        ctx.setLineDash([]);
        dBLabelData.map(x => {
            ctx.globalAlpha = 0.5;
            const label = `${x}dB`,
                posY = map(ascale(10 ** (x / 20)), 0, 1, canvas.height, 0);
            ctx.beginPath();
            ctx.lineTo(0, posY);
            ctx.lineTo(canvas.width, posY);
            ctx.stroke();
            ctx.globalAlpha = 1;
            ctx.textAlign = visualizerSettings.mirrorLabels ? 'end' : 'start'
            ctx.fillText(label, canvas.width * visualizerSettings.mirrorLabels, posY);
        });
        ctx.setLineDash([]);
        ctx.globalAlpha = 1;
        ctx.textAlign = 'start'
    }
    requestAnimationFrame(visualize);
    currentSampleRate = audioCtx.sampleRate;
}
// and here;s the additional functions that we can need for this visualization
function applyWindow(posX, windowType = 'Hann', windowParameter = 1, truncate = true, windowSkew = 0) {
    let x = windowSkew > 0 ? ((posX / 2 - 0.5) / (1 - (posX / 2 - 0.5) * 10 * (windowSkew ** 2))) / (1 / (1 + 10 * (windowSkew ** 2))) * 2 + 1 :
        ((posX / 2 + 0.5) / (1 + (posX / 2 + 0.5) * 10 * (windowSkew ** 2))) / (1 / (1 + 10 * (windowSkew ** 2))) * 2 - 1;

    if (truncate && Math.abs(x) > 1)
        return 0;

    switch (windowType.toLowerCase()) {
        default:
            return 1;
        case 'hanning':
        case 'cosine squared':
        case 'hann':
            return Math.cos(x * Math.PI / 2) ** 2;
        case 'raised cosine':
        case 'hamming':
            return 0.54 + 0.46 * Math.cos(x * Math.PI);
        case 'power of sine':
            return Math.cos(x * Math.PI / 2) ** windowParameter;
        case 'circle':
        case 'power of circle':
            return Math.sqrt(1 - (x ** 2)) ** windowParameter;
        case 'tapered cosine':
        case 'tukey':
            return Math.abs(x) <= 1 - windowParameter ? 1 :
                (x > 0 ?
                    (-Math.sin((x - 1) * Math.PI / windowParameter / 2)) ** 2 :
                    Math.sin((x + 1) * Math.PI / windowParameter / 2) ** 2);
        case 'blackman':
            return 0.42 + 0.5 * Math.cos(x * Math.PI) + 0.08 * Math.cos(x * Math.PI * 2);
        case 'nuttall':
            return 0.355768 + 0.487396 * Math.cos(x * Math.PI) + 0.144232 * Math.cos(2 * x * Math.PI) + 0.012604 * Math.cos(3 * x * Math.PI);
        case 'flat top':
        case 'flattop':
            return 0.21557895 + 0.41663158 * Math.cos(x * Math.PI) + 0.277263158 * Math.cos(2 * x * Math.PI) + 0.083578947 * Math.cos(3 * x * Math.PI) + 0.006947368 * Math.cos(4 * x * Math.PI);
        case 'kaiser':
            return Math.cosh(Math.sqrt(1 - (x ** 2)) * (windowParameter ** 2)) / Math.cosh(windowParameter ** 2);
        case 'gauss':
        case 'gaussian':
            return Math.exp(-(windowParameter ** 2) * (x ** 2));
        case 'cosh':
        case 'hyperbolic cosine':
            return Math.E ** (-(windowParameter ** 2) * (Math.cosh(x) - 1));
        case 'bartlett':
        case 'triangle':
        case 'triangular':
            return 1 - Math.abs(x);
        case 'poisson':
        case 'exponential':
            return Math.exp(-Math.abs(x * (windowParameter ** 2)));
        case 'hyperbolic secant':
        case 'sech':
            return 1 / Math.cosh(x * (windowParameter ** 2));
        case 'quadratic spline':
            return Math.abs(x) <= 0.5 ? -((x * Math.sqrt(2)) ** 2) + 1 : (Math.abs(x * Math.sqrt(2)) - Math.sqrt(2)) ** 2;
        case 'parzen':
            return Math.abs(x) > 0.5 ? -2 * ((-1 + Math.abs(x)) ** 3) : 1 - 24 * (Math.abs(x / 2) ** 2) + 48 * (Math.abs(x / 2) ** 3);
        case 'welch':
            return 1 - (x ** 2);
        case 'ogg':
        case 'vorbis':
            return Math.sin(Math.PI / 2 * Math.cos(x * Math.PI / 2) ** 2);
        case 'cascaded sine':
        case 'cascaded cosine':
        case 'cascaded sin':
        case 'cascaded cos':
            return 1 - Math.sin(Math.PI / 2 * Math.sin(x * Math.PI / 2) ** 2);
    }
}

function fscale(x, freqScale = 'logarithmic', freqSkew = 0.5) {
    switch (freqScale.toLowerCase()) {
        default:
            return x;
        case 'log':
        case 'logarithmic':
            return Math.log2(x);
        case 'mel':
            return Math.log2(1 + x / 700);
        case 'critical bands':
        case 'bark':
            return (26.81 * x) / (1960 + x) - 0.53;
        case 'equivalent rectangular bandwidth':
        case 'erb':
            return Math.log2(1 + 0.00437 * x);
        case 'cam':
        case 'cams':
            return Math.log2((x / 1000 + 0.312) / (x / 1000 + 14.675));
        case 'sinh':
        case 'arcsinh':
        case 'asinh':
            return Math.asinh(x / (10 ** (freqSkew * 4)));
        case 'shifted log':
        case 'shifted logarithmic':
            return Math.log2((10 ** (freqSkew * 4)) + x);
        case 'nth root':
            return x ** (1 / (11 - freqSkew * 10));
        case 'negative exponential':
            return -(2 ** (-x / (2 ** (7 + freqSkew * 8))));
        case 'adjustable bark':
            return (26.81 * x) / ((10 ** (freqSkew * 4)) + x);
        case 'period':
            return 1 / x;
    }
}

function invFscale(x, freqScale = 'logarithmic', freqSkew = 0.5) {
    switch (freqScale.toLowerCase()) {
        default:
            return x;
        case 'log':
        case 'logarithmic':
            return 2 ** x;
        case 'mel':
            return 700 * ((2 ** x) - 1);
        case 'critical bands':
        case 'bark':
            return 1960 / (26.81 / (x + 0.53) - 1);
        case 'equivalent rectangular bandwidth':
        case 'erb':
            return (1 / 0.00437) * ((2 ** x) - 1);
        case 'cam':
        case 'cams':
            return (14.675 * (2 ** x) - 0.312) / (1 - (2 ** x)) * 1000;
        case 'sinh':
        case 'arcsinh':
        case 'asinh':
            return Math.sinh(x) * (10 ** (freqSkew * 4));
        case 'shifted log':
        case 'shifted logarithmic':
            return (2 ** x) - (10 ** (freqSkew * 4));
        case 'nth root':
            return x ** ((11 - freqSkew * 10));
        case 'negative exponential':
            return -Math.log2(-x) * (2 ** (7 + freqSkew * 8));
        case 'adjustable bark':
            return (10 ** (freqSkew * 4)) / (26.81 / x - 1);
        case 'period':
            return 1 / x;
    }
}

function generateFreqBands(N = 128, low = 20, high = 20000, freqScale, freqSkew, bandwidth = 0.5) {
    let freqArray = [];
    for (let i = 0; i < N; i++) {
        freqArray.push({
            lo: invFscale(map(i - bandwidth, 0, N - 1, fscale(low, freqScale, freqSkew), fscale(high, freqScale, freqSkew)), freqScale, freqSkew),
            ctr: invFscale(map(i, 0, N - 1, fscale(low, freqScale, freqSkew), fscale(high, freqScale, freqSkew)), freqScale, freqSkew),
            hi: invFscale(map(i + bandwidth, 0, N - 1, fscale(low, freqScale, freqSkew), fscale(high, freqScale, freqSkew)), freqScale, freqSkew)
        });
    }
    return freqArray;
}

function generateOctaveBands(bandsPerOctave = 12, lowerNote = 4, higherNote = 123, detune = 0, tuningFreq = 440, bandwidth = 0.5) {
    const tuningNote = isFinite(Math.log2(tuningFreq)) ? Math.round((Math.log2(tuningFreq) - 4) * 12) * 2 : 0,
        root24 = 2 ** (1 / 24),
        c0 = tuningFreq * root24 ** -tuningNote, // ~16.35 Hz
        groupNotes = 24 / bandsPerOctave;
    let bands = [];
    for (let i = Math.round(lowerNote * 2 / groupNotes); i <= Math.round(higherNote * 2 / groupNotes); i++) {
        bands.push({
            lo: c0 * root24 ** ((i - bandwidth) * groupNotes + detune),
            ctr: c0 * root24 ** (i * groupNotes + detune),
            hi: c0 * root24 ** ((i + bandwidth) * groupNotes + detune)
        });
    }
    return bands;
}

function ascale(x) {
    if (visualizerSettings.useDecibels)
        return map(20 * Math.log10(x), visualizerSettings.minDecibels, visualizerSettings.maxDecibels, 0, 1);
    else
        return map(x ** (1 / visualizerSettings.gamma), !visualizerSettings.useAbsolute * (10 ** (visualizerSettings.minDecibels / 20)) ** (1 / visualizerSettings.gamma), (10 ** (visualizerSettings.maxDecibels / 20)) ** (1 / visualizerSettings.gamma), 0, 1);
}

function parseList(string) {
    return string.split(',').map(x => isNaN(x) ? 0 : parseFloat(x));
}

function updateSpectrumVisualization(data, inAudioContext = false) {
    if (currentSpectrum.length !== data.length)
        currentSpectrum.length = data.length;
    if (currentSpectrum.length !== peaks.length || currentSpectrum.length !== peakHolds.length) {
        peaks.length = currentSpectrum.length;
        peakHolds.length = currentSpectrum.length;
    }
    const factor = inAudioContext ? 60 / audioCtx.sampleRate : 1,
        holdFactor = inAudioContext ? audioCtx.sampleRate / 60 : 1,
        smoothingTimeConstant = (visualizerSettings.smoothingTimeConstant / 100) ** factor,
        peakDecayTimeConstant = (visualizerSettings.peakDecay / 100) ** factor;
    for (let i = 0; i < data.length; i++) {
        currentSpectrum[i] = isFinite(currentSpectrum[i]) ? visualizerSettings.useAverageSmoothing ? data[i] * (1 - smoothingTimeConstant) + currentSpectrum[i] * smoothingTimeConstant : Math.max(data[i], currentSpectrum[i] * smoothingTimeConstant) : data[i];
        const peakValue = visualizerSettings.useActualPeak ? data[i] : currentSpectrum[i];
        if (peakValue >= peaks[i] || !isFinite(peaks[i])) {
            peaks[i] = peakValue;
            peakHolds[i] = visualizerSettings.peakHold * holdFactor;
        } else if (peakHolds[i] > 0)
            peakHolds[i] = Math.min(peakHolds[i] - 1, visualizerSettings.peakHold * holdFactor);
        else
            peaks[i] *= peakDecayTimeConstant;
    }
}
```
