// audio-worklet.js
// // @ts-ignore
// class GainUpProcessor extends AudioWorkletProcessor {
//   constructor(context) {
//     super(context, 'GainUpProcessor')
//   }
//   process(inputs, outputs, parameters) {
//     const input = inputs[0]
//     const output = outputs[0]

//     const gainFactor = 10 // 10 dB gain

//     for (let channel = 0; channel < input.length; channel++) {
//       for (let i = 0; i < input[channel].length; i++) {
//         output[channel][i] = input[channel][i] * gainFactor
//       }
//     }

//     return true
//   }
// }

// This is "processor.js" file, evaluated in AudioWorkletGlobalScope upon
// audioWorklet.addModule() call in the main global scope.
class MyWorkletProcessor extends AudioWorkletProcessor {
  constructor() {
    super()
  }

  process(inputs, outputs, parameters) {
    const input = inputs[0]
    const output = outputs[0]

    const gainFactor = 10 // 10 dB gain

    for (let channel = 0; channel < input.length; channel++) {
      for (let i = 0; i < input[channel].length; i++) {
        output[channel][i] = input[channel][i] * gainFactor
      }
    }

    return true
  }
}

registerProcessor('my-worklet-processor', MyWorkletProcessor)
