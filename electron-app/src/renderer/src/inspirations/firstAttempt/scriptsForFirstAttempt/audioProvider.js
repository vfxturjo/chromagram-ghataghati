{
  /* <script id="AudioProvider" type="worklet"> */
}
class AudioProvider extends AudioWorkletProcessor {
  constructor() {
    super()
    this.dataArrays = []
    this.bufferSize = 32768 // can handle more than 32768 samples of PCM data unlike in AnalyserNode.getFloatTimeDomainData, which is capped at 32768 samples
    this.bufferIdx = 0
    this.currentTimeInSamples = 0
    this.port.onmessage = (e) => {
      const audioChunks = [],
        retrievalWindowSize = e.data
          ? Math.min(this.bufferSize, currentFrame - this.currentTimeInSamples)
          : this.bufferSize,
        timeOffset = this.bufferSize - retrievalWindowSize
      for (let channelIdx = 0; channelIdx < this.dataArrays.length; channelIdx++) {
        audioChunks[channelIdx] = []
        for (let i = 0; i < this.dataArrays[channelIdx].length - timeOffset; i++) {
          const data =
            this.dataArrays[channelIdx][
              (((this.bufferIdx + i + timeOffset) % this.bufferSize) + this.bufferSize) %
                this.bufferSize
            ]
          audioChunks[channelIdx][i] = data !== undefined ? data : 0
        }
      }
      this.port.postMessage({ currentChunk: audioChunks })
      this.currentTimeInSamples = currentFrame
    }
  }

  process(inputs, _, _2) {
    if (inputs[0].length <= 0) return true
    this.dataArrays.length = inputs[0].length
    for (let i = 0; i < this.dataArrays.length; i++) {
      if (this.dataArrays[i] === undefined) this.dataArrays[i] = new Array(this.bufferSize)
      else {
        this.dataArrays[i].length = this.bufferSize
      }
    }

    for (let i = 0; i < inputs[0][0].length; i++) {
      this.bufferIdx = Math.min(this.bufferIdx, this.bufferSize - 1)
      for (let channelIdx = 0; channelIdx < inputs[0].length; channelIdx++) {
        this.dataArrays[channelIdx][this.bufferIdx] = inputs[0][channelIdx][i]
      }
      this.bufferIdx =
        (((this.bufferIdx + 1) % this.bufferSize) + this.bufferSize) % this.bufferSize
    }
    return true
  }
}

registerProcessor('audio-provider', AudioProvider)
{
  /* </script> */
}
