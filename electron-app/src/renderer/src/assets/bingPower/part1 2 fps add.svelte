<!-- Chromagram.svelte -->
<script>
  import { onMount } from 'svelte'

  let audioContext
  let analyser
  let chromagram = new Array(12).fill(0)
  let fps = 30 // Default frame rate

  async function startAudio() {
    audioContext = new AudioContext()
    const stream = await navigator.mediaDevices.getUserMedia({ audio: true })
    const source = audioContext.createMediaStreamSource(stream)
    analyser = audioContext.createAnalyser()
    source.connect(analyser)
    analyser.fftSize = 4096 // Adjust as needed
  }

  function calculateChromagram() {
    const dataArray = new Uint8Array(analyser.frequencyBinCount)
    analyser.getByteFrequencyData(dataArray)

    // Map frequency data to pitch classes
    dataArray.forEach((value, bin) => {
      const frequency = (bin * audioContext.sampleRate) / analyser.fftSize
      const pitchClassIndex = Math.floor((12 * Math.log2(frequency / 440)) % 12)
      chromagram[pitchClassIndex] += value
    })

    // Normalize the chromagram (optional)
    const totalEnergy = chromagram.reduce((sum, val) => sum + val, 0)
    chromagram = chromagram.map((val) => val / totalEnergy)

    // Call the visualization function
    visualizeChromagram(chromagram)
  }

  // Function to update chromagram at specified frame rate
  let interval
  function updateChromagram() {
    clearInterval(interval)
    interval = setInterval(calculateChromagram, 1000 / fps)
  }

  // // // // // // // //
  // // // //
  //#region Drawing the chromagram
  let canvasChroma12
  let ctx
  let canvasWidth
  let canvasHeight
  let barWidth

  onMount(() => {
    ctx = canvasChroma12.getContext('2d')
    canvasWidth = canvasChroma12.width
    canvasHeight = canvasChroma12.height
    // barWidth = canvasWidth / chromagramData.length
    barWidth = canvasWidth / 12
    console.log(ctx, canvasHeight, canvasWidth, barWidth)
  })

  // Visualization function
  function visualizeChromagram(chromagramData) {
    // Clear the canvas
    ctx.clearRect(0, 0, canvasWidth, canvasHeight)
    // Draw chromagram bars
    chromagramData.forEach((value, index) => {
      const barHeight = value * canvasHeight
      ctx.fillStyle = 'red'

      ctx.fillRect(index * barWidth, canvasHeight - barHeight, barWidth, barHeight)
    })
  }
  //#endregion
</script>

<button on:click={startAudio}>Start Audio</button>
<button on:click={calculateChromagram}>Calculate Chromagram</button>

<label>
  Frame Rate (fps):
  <input type="range" min="1" max="60" step="1" bind:value={fps} on:click={updateChromagram} />
  {fps} fps
</label>

<canvas
  bind:this={canvasChroma12}
  id="chromagram-canvas"
  width="480"
  height="320"
  style="border: 2px solid white;"
></canvas>

<!-- Display the chromagram values -->
{#each chromagram as value, index}
  <p>{index}: {value.toFixed(4)}</p>
{/each}
