// AudioInput.js
export class AudioInput {
  constructor() {
    // @ts-ignore
    this.audioContext = new (window.AudioContext || window.webkitAudioContext)()
    this.microphone = null
    this.scriptProcessor = null
    this.isCapturing = false
  }

  async start() {
    if (this.isCapturing) {
      console.warn('Audio capture is already running.')
      return
    }

    try {
      const stream = await navigator.mediaDevices.getUserMedia({ audio: true })
      this.microphone = this.audioContext.createMediaStreamSource(stream)
      this.scriptProcessor = this.audioContext.createScriptProcessor(4096, 1, 1)

      this.microphone.connect(this.scriptProcessor)
      this.scriptProcessor.connect(this.audioContext.destination)
      this.scriptProcessor.onaudioprocess = this.processAudio.bind(this)

      this.isCapturing = true
    } catch (error) {
      console.error('Error starting audio capture:', error)
    }
  }

  stop() {
    if (!this.isCapturing) {
      console.warn('Audio capture is not running.')
      return
    }

    if (this.scriptProcessor) {
      this.scriptProcessor.disconnect()
      this.scriptProcessor.onaudioprocess = null
      this.scriptProcessor = null
    }

    if (this.microphone) {
      this.microphone.disconnect()
      this.microphone = null
    }

    this.audioContext.close()
    this.isCapturing = false
  }

  processAudio(event) {
    const inputBuffer = event.inputBuffer
    const inputData = inputBuffer.getChannelData(0)
    // Emit the raw audio data for further processing
    this.emit('audio-data', inputData)
  }

  on(event, handler) {
    this.events = this.events || {}
    this.events[event] = this.events[event] || []
    this.events[event].push(handler)
  }

  emit(event, data) {
    if (this.events && this.events[event]) {
      this.events[event].forEach((handler) => handler(data))
    }
  }
}
