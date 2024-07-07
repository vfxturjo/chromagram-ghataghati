<script lang="ts">
  let counter = 0
</script>

This is a counter.
<button type="button" class="btn variant-filled" on:click={() => counter++}
  >Clicked: {counter}</button
>
