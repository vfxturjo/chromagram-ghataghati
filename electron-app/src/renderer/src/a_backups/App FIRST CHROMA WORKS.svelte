<script lang="ts">
  import { onMount } from 'svelte'
  // import { Chromagram } from './Chromagram' // Assuming Chromagram.ts is in the same directory
  import { Chromagram } from '../assets/AKKORDER/Chromagram'
  import Range from '../components/Slider.svelte'

  let audioContext
  let microphone
  let analyser
  let dataArray
  let canvas, ctx
  let chromagram
  let frameSize = 4096 // Example frame size for Chromagram
  let samplingFrequency = 44100 // Example sampling frequency
  let isAudioSetup = false
  let power = 1

  async function setupAudio() {
    audioContext = new (window.AudioContext || window.webkitAudioContext)()
    const stream = await navigator.mediaDevices.getUserMedia({ audio: true })
    microphone = audioContext.createMediaStreamSource(stream)
    analyser = audioContext.createAnalyser()
    analyser.fftSize = frameSize * 2
    dataArray = new Float32Array(analyser.fftSize)
    microphone.connect(analyser)

    chromagram = new Chromagram(frameSize, samplingFrequency)
    isAudioSetup = true
    draw()
  }

  function draw() {
    if (!isAudioSetup) return
    requestAnimationFrame(draw)

    analyser.getFloatTimeDomainData(dataArray)
    chromagram.processAudioFrame(dataArray)

    if (chromagram.isReady()) {
      let chromaData = chromagram.getChromagram()

      // MY THINGGGGG

      chromaData = normalizeData(chromaData, power)

      // MY THINGGGGG

      drawChromagram(chromaData)
    }
  }

  function drawChromagram(chromaData) {
    const width = canvas.width
    const height = canvas.height
    const barWidth = width / chromaData.length

    ctx.clearRect(0, 0, width, height)

    chromaData.forEach((value, index) => {
      const barHeight = value * height
      ctx.fillStyle = `hsl(${(index / chromaData.length) * 360}, 100%, 50%)`
      ctx.fillRect(index * barWidth, height - barHeight, barWidth, barHeight)
    })
  }

  onMount(() => {
    canvas = document.querySelector('#chromagramCanvas')
    ctx = canvas.getContext('2d')
    setupAudio()
  })

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
<canvas id="chromagramCanvas" width="640" height="360"></canvas>
<div class:purple-theme={theme === 'purple'} style="width: 100%;">
  <label for="basic-range" style="width: 100%;">Range Label</label>
  <Range on:change={(e) => (power = e.detail.value)} id="basic-slider" />
</div>
