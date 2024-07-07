import { join } from 'path'
import type { Config } from 'tailwindcss'
import forms from '@tailwindcss/forms'
import typography from '@tailwindcss/typography'

// 1. Import the Skeleton plugin
import { skeleton } from '@skeletonlabs/tw-plugin'

const config = {
  darkMode: 'class',
  content: [
    './src/**/*.{html,js,svelte,ts}',
    // join(require.resolve('@skeletonlabs/skeleton'), '../**/*.{html,js,svelte,ts}')
    // join('./node_modules/@skeletonlabs/skeleton', '../**/*.{html,js,svelte,ts}')
    './node_modules/@skeletonlabs/**/*.{html,js,svelte,ts}'
  ],

  theme: {
    extend: {}
  },

  plugins: [
    forms,
    typography,
    skeleton({
      themes: {
        preset: [
          {
            name: 'skeleton',
            enhancements: true
          }
        ]
      }
    })
  ]
} satisfies Config

export default config
