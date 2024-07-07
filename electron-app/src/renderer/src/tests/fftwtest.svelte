<script lang="ts">
  import { onMount } from 'svelte'
  // import * as FFTWModule from '../_jsAssets/FFTW/FFTW'
  // import '../_jsAssets/FFTW/FFTW.js'
  // import '../_jsAssets/FFTW/FFT.js'
  import { FFTW } from '../_jsAssets/FFTW/FFT'

  let audioContext
  let fftSize = 2048
  let canvas, canvasCtx, canvas2, canvas2Ctx
  let FFTWworker = new FFTW(fftSize)

  // AudioWorkletProcessor class should be defined in a separate file 'audio-processor.js'
  const audioWorkletProcessorUrl = 'audio-processor.js'

  onMount(async () => {
    canvas = document.getElementById('fftCanvas')
    canvasCtx = canvas.getContext('2d')
    canvas2 = document.getElementById('byteFrequencyCanvas')
    canvas2Ctx = canvas2.getContext('2d')

    if (!navigator.mediaDevices.getUserMedia) {
      console.error('getUserMedia not supported on your browser.')
      return
    }

    try {
      audioContext = new AudioContext()
      await audioContext.audioWorklet.addModule(audioWorkletProcessorUrl)
      const stream = await navigator.mediaDevices.getUserMedia({ audio: true })
      const source = audioContext.createMediaStreamSource(stream)
      const workletNode = new AudioWorkletNode(audioContext, 'audio-processor')

      source.connect(workletNode)
      workletNode.connect(audioContext.destination)

      workletNode.port.onmessage = (event) => {
        // Assuming the processor sends audio data in event.data
        const inputData = event.data

        // Calculate FFT using fftw-js
        let fftOutput = FFTWworker.forward(inputData)
        visualizeFFT(fftOutput, canvasCtx)

        // Calculate FFT using getByteFrequencyData
        let analyser = audioContext.createAnalyser()
        source.connect(analyser)
        let dataArray = new Uint8Array(analyser.frequencyBinCount)
        analyser.getByteFrequencyData(dataArray)
        visualizeByteFrequency(dataArray, canvas2Ctx)
      }
    } catch (error) {
      console.error('The following error occurred: ' + error)
    }
  })

  function visualizeFFT(fftData, ctx) {
    ctx.clearRect(0, 0, canvas.width, canvas.height)
    ctx.fillStyle = 'rgb(0, 0, 0)'
    ctx.fillRect(0, 0, canvas.width, canvas.height)
    ctx.lineWidth = 2
    ctx.strokeStyle = 'rgb(0, 255, 0)'
    ctx.beginPath()

    let sliceWidth = (canvas.width * 1.0) / fftData.length
    let x = 0

    for (let i = 0; i < fftData.length; i++) {
      let v = fftData[i] / 128.0
      let y = (v * canvas.height) / 2

      if (i === 0) {
        ctx.moveTo(x, y)
      } else {
        ctx.lineTo(x, y)
      }

      x += sliceWidth
    }

    ctx.lineTo(canvas.width, canvas.height / 2)
    ctx.stroke()
  }

  function visualizeByteFrequency(dataArray, ctx) {
    ctx.clearRect(0, 0, canvas2.width, canvas2.height)
    ctx.fillStyle = 'rgb(0, 0, 0)'
    ctx.fillRect(0, 0, canvas2.width, canvas2.height)
    let barWidth = (canvas2.width / dataArray.length) * 2.5
    let barHeight
    let x = 0

    for (let i = 0; i < dataArray.length; i++) {
      barHeight = dataArray[i]

      ctx.fillStyle = 'rgb(' + (barHeight + 100) + ',50,50)'
      ctx.fillRect(x, canvas2.height - barHeight / 2, barWidth, barHeight / 2)

      x += barWidth + 1
    }
  }
</script>

<!-- <svelte:head>
  <script
    src="../_jsAssets/FFTW/FFTW.js"
    on:load={() => console.log('script:load FFTW main')}
  ></script>
  <script
    src="../_jsAssets/FFTW/FFT.js"
    on:load={() => console.log('script:load FFTW wrapper')}
  ></script>
</svelte:head> -->

<canvas id="fftCanvas" width="640" height="480"></canvas>
<canvas id="byteFrequencyCanvas" width="640" height="480"></canvas>
