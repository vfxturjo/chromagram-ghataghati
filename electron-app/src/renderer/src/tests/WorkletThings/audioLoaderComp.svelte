<script lang="ts">
  export let stream
  var processorPath = '/src/assets/Worklet-Processors/processor.js'
  var testing_situation
  var context
  var resulting_thing

  async function startMicrophone() {
    // Prompt the user to use their microphone.
    stream = await navigator.mediaDevices.getUserMedia({
      audio: true
    })

    context = new AudioContext()
    var source = context.createMediaStreamSource(stream)

    // Load and execute the module script.
    await context.audioWorklet.addModule(processorPath)
    // Create an AudioWorkletNode. The name of the processor is the
    // one passed to registerProcessor() in the module script.
    const processor = new AudioWorkletNode(context, 'processor')

    source.connect(processor).connect(resulting_thing)
    console.log('Your microphone audio is being used.')
  }

  function stopMicrophone() {
    // Stop the stream.
    stream.getTracks().forEach((track) => track.stop())
    context.close()
    console.log('Your microphone audio is not used anymore.')
  }
</script>

<button on:click={startMicrophone}>Start microphone</button>
<button on:click={stopMicrophone}>Stop microphone</button>
<button
  on:click={() => {
    console.log(resulting_thing)
  }}>LOG STREAM thingg</button
>
