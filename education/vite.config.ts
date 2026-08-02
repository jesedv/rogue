import { defineConfig } from "vite";
import wasm from "vite-plugin-wasm";
import topLevelAwait from "vite-plugin-top-level-await";

export default defineConfig({
  plugins: [wasm(), topLevelAwait()],
  server: {
    port: 5173,
    fs: { allow: [".."] },
    configureServer(server) {
      server.middlewares.use((req, _res, next) => {
        if (req.url! === "/education" || req.url!.startsWith("/education/")) {
          req.url = "/education.html";
        } else if (req.url! === "/production" || req.url!.startsWith("/production/")) {
          req.url = "/production.html";
        } else if (req.url! === "/") {
          req.url = "/index.html";
        }
        next();
      });
    },
  },
  optimizeDeps: {
    exclude: ["rogue-wasm"],
  },
  build: {
    rollupOptions: {
      input: {
        index: new URL("./index.html", import.meta.url).pathname,
        education: new URL("./education.html", import.meta.url).pathname,
        production: new URL("./production.html", import.meta.url).pathname,
      },
    },
  },
});
