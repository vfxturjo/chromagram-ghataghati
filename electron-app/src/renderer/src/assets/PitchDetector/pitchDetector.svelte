<script>
  import { onMount } from 'svelte'
  import Pitchfinder from 'pitchfinder'

  let audioContext
  let micStream
  let audioStreamSource
  let pitchDetectorYIN
  let pitchDetectorAMDF
  let animationFrameId
  let pitch1
  let pitch2
  let dataArray
  let analyserNode

  onMount(async () => {
    try {
      // Initialize audio context and pitch detection
      audioContext = new AudioContext()
      pitchDetectorYIN = Pitchfinder.YIN()
      pitchDetectorAMDF = Pitchfinder.AMDF()

      // Request access to the microphone
      const stream = await navigator.mediaDevices.getUserMedia({ audio: true })
      micStream = stream
      audioStreamSource = audioContext.createMediaStreamSource(stream)

      // Create an AnalyserNode to process audio data
      analyserNode = audioContext.createAnalyser()
      analyserNode.fftSize = 2048
      audioStreamSource.connect(analyserNode)

      // Process audio data
      const bufferLength = analyserNode.frequencyBinCount
      dataArray = new Float32Array(bufferLength)

      // Start the animation frame loop
      animationFrameId = requestAnimationFrame(updatePitch)
    } catch (error) {
      console.error('Error accessing the microphone', error)
    }
  })

  function updatePitch() {
    analyserNode.getFloatTimeDomainData(dataArray)
    pitch1 = pitchDetectorYIN(dataArray)
    pitch2 = pitchDetectorAMDF(dataArray)

    // Continue the animation frame loop
    animationFrameId = requestAnimationFrame(updatePitch)
  }
</script>

<p>Pitch1 is: {Math.round(pitch1)}</p>
<p>Pitch2 is: {Math.round(pitch2)}</p>
