import { defineConfig } from "vite";
import wasm from "vite-plugin-wasm";
import topLevelAwait from "vite-plugin-top-level-await";

export default defineConfig({
  plugins: [wasm(), topLevelAwait()],
  server: {
    port: 5173,
    fs: { allow: [".."] },
  },
  optimizeDeps: {
    exclude: ["rogue-wasm"],
  },
  build: {
    rollupOptions: {
      input: {
        index: new URL("./index.html", import.meta.url).pathname,
        production: new URL("./production.html", import.meta.url).pathname,
      },
    },
  },
});
