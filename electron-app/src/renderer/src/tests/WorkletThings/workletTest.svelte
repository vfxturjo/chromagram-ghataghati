<script lang="ts">
  let stream
  import AudioLoaderComp from './audioLoaderComp.svelte'
</script>

<AudioLoaderComp bind:stream></AudioLoaderComp>

<button
  on:click={() => {
    console.log(stream)
  }}>LOG STREAM</button
>
