<script lang="ts">
  import { CircularArray } from '../assets/ringBuffer/ringBuffer'

  var justNumberArray = [1, 2, 3, 4, 5]
  var testArray = new CircularArray(10, ...justNumberArray)
  console.log(testArray)

  function addElementsInTestArray() {
    justNumberArray = justNumberArray.map((i) => {
      return i + 1
    })

    testArray.push(...justNumberArray)
    console.log(testArray)
    testArray = testArray
  }
</script>

<button on:click={addElementsInTestArray}>CHANGE</button>

{#each testArray as item}
  {item},
{:else}
  <!-- empty list -->
{/each}
