<script lang="ts">
  import { onMount } from 'svelte'
  // import FFTWModule from '../_jsAssets/FFTW/FFTW'
  // import FFT from 'fftw-js/src/main'
  import { Chromagram } from '../assets/AKKORDER/Chromagram mod'
  import Range from '../components/Slider.svelte'

  let audioContext
  let microphone
  let analyser
  let dataArray
  let canvas, ctx
  let waterfallCanvas, waterfallCtx
  let chromagram
  let frameSize = 4096 // Example frame size for Chromagram
  let samplingFrequency = 44100 // Example sampling frequency
  let isAudioSetup = false
  let power = 8
  let historyData = []
  let chromaData
  let stream

  let mediaRecorder
  let chunks = []

  onMount(() => {
    ctx = canvas.getContext('2d')
    waterfallCtx = waterfallCanvas.getContext('2d')
  })

  async function setupAudio() {
    audioContext = new (window.AudioContext || window.webkitAudioContext)()
    stream = await navigator.mediaDevices.getUserMedia({ audio: true })
    microphone = audioContext.createMediaStreamSource(stream)
    analyser = audioContext.createAnalyser()
    analyser.fftSize = frameSize * 2
    dataArray = new Float32Array(analyser.fftSize)
    microphone.connect(analyser)

    mediaRecorder = new MediaRecorder(stream)
    mediaRecorder.onDataAvailable = (event) => {
      console.log(event)

      chunks.push(event.data)
    }

    chromagram = new Chromagram(frameSize, samplingFrequency, 100)
    isAudioSetup = true
    draw()
  }

  function draw() {
    if (!isAudioSetup) return
    requestAnimationFrame(draw)

    analyser.getFloatTimeDomainData(dataArray)

    chromagram.processAudioFrame(dataArray)

    if (chromagram.isReady()) {
      chromaData = chromagram.getChromagram()
      chromaData = normalizeData(chromaData, power)
      drawChromagram()
      // drawWaterfallChroma()
    }
  }

  function drawChromagram() {
    const width = canvas.width
    const height = canvas.height
    const totalBars = chromaData.length // Now we have 120 bars
    const barWidth = width / totalBars

    ctx.clearRect(0, 0, width, height)

    chromaData.forEach((value, index) => {
      const barHeight = value * height
      // We can use a more detailed color scale to represent the 120 bars
      ctx.fillStyle = `hsl(${(index / totalBars) * 360}, 100%, 50%)`
      ctx.fillRect(index * barWidth, height - barHeight, barWidth, barHeight)
    })
  }

  function drawWaterfallChroma() {
    const width = waterfallCanvas.width
    const height = waterfallCanvas.height
    const totalBars = chromaData.length
    const barHeight = height / totalBars

    // Shift history data to the left
    if (historyData.length > 0) {
      waterfallCtx.putImageData(historyData[0], -1, 0)
      historyData[0] = waterfallCtx.getImageData(0, 0, width, height)
    }

    // Clear the rightmost column to draw the new data
    waterfallCtx.clearRect(width - 1, 0, 1, height)

    // Draw the current chromaData on the rightmost part of the canvas
    chromaData.forEach((value, index) => {
      // Use the index to set the y position of the bar
      const yPos = height - index * barHeight
      // Use the value to set the brightness of the rectangle's color
      const brightness = value * 100 // Assuming 'value' is normalized between 0 and 1
      waterfallCtx.fillStyle = `hsl(0, 100%, ${brightness}%)`
      waterfallCtx.fillRect(width - 1, yPos, 1, barHeight)
    })

    // Save the current frame to history
    historyData.unshift(waterfallCtx.getImageData(0, 0, width, height))
  }

  function normalizeData(data, power) {
    // Raise each element to the specified power and sum them up
    const total = data.reduce((sum, value) => sum + Math.pow(value, power), 0)

    data = data.map((value) => Math.pow(value, power) / total)
    const maxValue = Math.max(...data)
    // return data.map((value) => value / maxValue)

    // Normalize each element by dividing by the total sum
    return data.map((value) => value / maxValue)
  }

  // SLIDER THINGS
  let theme = 'default'
</script>

<button on:click={setupAudio}>Start Microphone</button>
<button on:click={() => console.log(chromagram.getChromagram())}>show chromaData</button>
<canvas id="chromagramCanvas" width="640" height="100" bind:this={canvas}></canvas>
<canvas id="waterfallCanvas" width="640" height="360" bind:this={waterfallCanvas}></canvas>
<div class:purple-theme={theme === 'purple'} style="width: 100%;">
  <label for="basic-range" style="width: 100%;">Range Label</label>
  <Range initialValue={power} on:change={(e) => (power = e.detail.value)} id="basic-slider" />
</div>
<button
  on:click={() => {
    // console.log(dataArray)
    console.log(chunks)
  }}>OUTPUT TIME DOMAIN DATA FROM BOTH WORLD</button
>
