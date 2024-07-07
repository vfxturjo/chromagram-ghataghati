<!--
Created: Sun May 26 2024 02:07:29 GMT+0600 (Bangladesh Standard Time)
Modified: Sun May 26 2024 03:02:16 GMT+0600 (Bangladesh Standard Time)
-->

# Drawing the chromagram data on canvas

```svelte
<script>
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

<!-- main Chromagram canvas -->
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

```

# Loading FFTW and RFFTW

```js
import FFT, RFFT from 'fftw-js/src/main'
```

# Chromagram Update Frequency (FPS)

```svelte
<script>
  // Function to update chromagram at specified frame rate
  let fps = 30 // Default frame rate
  let interval
  function updateChromagram() {
    clearInterval(interval)
    interval = setInterval(calculateChromagram, 1000 / fps)
  }
</script>

<label>
  Frame Rate (fps):
  <input type="range" min="1" max="60" step="1" bind:value={fps} on:click={updateChromagram} />
  {fps} fps
</label>
```
