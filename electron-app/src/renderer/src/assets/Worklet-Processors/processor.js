// processor.js
// This file is evaluated in the audio rendering thread
// upon context.audioWorklet.addModule() call.

export var testingVar = 0
var remaining = 100

class Processor extends AudioWorkletProcessor {
  process([input], [output]) {
    // Copy inputs to outputs.
    console.log(input)
    // output[0].set(input[0])

    // testingVar = input[0]
    // output[0] = [0, 0, 0, 0, 0, 0]
    remaining--
    if (remaining < 1) {
      return false
    } else {
      return true
    }
    // return true
  }
}

registerProcessor('processor', Processor)
