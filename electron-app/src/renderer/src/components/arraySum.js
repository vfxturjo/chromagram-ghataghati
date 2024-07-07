export function arraySum(arr) {
  let sum = 0

  function sumElements(element) {
    if (Array.isArray(element)) {
      element.forEach(sumElements)
    } else {
      sum += element
    }
  }

  arr.forEach(sumElements)
  return sum
}
