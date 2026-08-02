import init, { NlsSim, version } from "../pkg/rogue_wasm.js";

const canvas = document.getElementById("c") as HTMLCanvasElement;
const ctx = canvas.getContext("2d", { alpha: false })!;
const statsEl = document.getElementById("stats")!;
const etaEl = document.getElementById("eta")!;
const scenarioSel = document.getElementById("scenario") as HTMLSelectElement;
const resSel = document.getElementById("res") as HTMLSelectElement;
const speedEl = document.getElementById("speed") as HTMLInputElement;
const speedVal = document.getElementById("speedVal")!;

let sim: NlsSim | null = null;
let paused = false;
let nx = 512;
let lx = 80.0;
let frames = 0;
let lastFps = 0;
let fpsT0 = performance.now();
let subSteps = 1;

const img = { data: null as ImageData | null, w: 0, h: 0 };
let yoff = 0;

function resize() {
  const dpr = Math.min(devicePixelRatio, 2);
  const rect = canvas.parentElement!.getBoundingClientRect();
  canvas.width = Math.floor(rect.width * dpr);
  canvas.height = Math.floor(rect.height * dpr);
}

function ensureImage(w: number, h: number) {
  if (!img.data || img.w !== w || img.h !== h) {
    img.data = ctx.createImageData(w, h);
    img.w = w;
    img.h = h;
  }
  return img.data;
}

function heat(amp: Float32Array, sigma: number) {
  const image = ensureImage(nx, nx);
  const d = image.data;
  const inv = sigma > 0 ? 1 / sigma : 1;
  for (let i = 0; i < nx * nx; i++) {
    const v = Math.min(amp[i] * inv * 0.72, 3.0);
    const p = i * 4;
    if (v > 1.0) {
      const t = (v - 1.0) / 2.0;
      d[p] = 255;
      d[p + 1] = Math.floor(255 * (1 - t * 0.85));
      d[p + 2] = Math.floor(255 * (1 - t * 0.55));
    } else {
      d[p] = Math.floor(16 + 120 * v);
      d[p + 1] = Math.floor(16 + 170 * v);
      d[p + 2] = Math.floor(40 + 200 * v);
    }
    d[p + 3] = 255;
  }
  return image;
}

function draw() {
  if (!sim) return;
  const amp = sim.amplitudes();
  const diag = sim.diagnostics() as { mass: number; h1: number; energy: number; sup: number };
  const sigma = Math.sqrt(diag.mass / nx); // rms |ψ| over the box
  const image = heat(amp, sigma);

  const off = document.createElement("canvas");
  off.width = nx;
  off.height = nx;
  const octx = off.getContext("2d")!;
  octx.putImageData(image, 0, 0);

  // scroll the field vertically to visualize the moving envelope
  yoff = (yoff + 1) % nx;
  ctx.fillStyle = "#000";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.imageSmoothingEnabled = true;
  ctx.save();
  ctx.translate(0, canvas.height);
  ctx.scale(canvas.width / nx, -canvas.height / nx);
  ctx.drawImage(off, 0, yoff - nx);
  ctx.drawImage(off, 0, yoff);
  ctx.restore();

  const blow = sim.blow_up_state() as {
    active: boolean;
    guaranteed_by_energy: boolean;
    growth_rate: number;
    eta: number | null;
    eta_r2: number;
    concentration: number;
    virial_sign: number;
    supercritical: boolean;
  };
  const rogue = sim.rogue_stats() as {
    sigma: number;
    crest_factor: number;
    event_count: number;
    kurtosis: number;
  };

  let etaTxt = "—";
  let cls = "";
  if (blow.eta != null) {
    etaTxt = `ETA ${blow.eta.toFixed(3)} s (r²=${blow.eta_r2.toFixed(2)})`;
    cls = blow.eta < 0.15 ? "danger" : "warn";
  }
  etaEl.className = cls;
  etaEl.textContent = etaTxt;

  const lines = [
    `rogue ${version()}`,
    `fps  ${lastFps.toFixed(0)} · t ${sim.time().toFixed(2)}`,
    `M    ${diag.mass.toFixed(4)}`,
    `E    ${diag.energy.toFixed(4)}`,
    `H¹   ${diag.h1.toFixed(4)}`,
    `crest  ${diag.sup.toFixed(3)} ×σ ${rogue.crest_factor.toFixed(2)}`,
    `rogue events  ${rogue.event_count}`,
    `kurtosis ${rogue.kurtosis.toFixed(2)} <span class="kurt">`,
    `blow-up ${blow.active ? "ACTIVE" : "off"} ${blow.guaranteed_by_energy ? "· E<0" : ""}`,
  ];
  lines[lines.length - 2] = `kurtosis ${rogue.kurtosis.toFixed(2)}`;
  statsEl.innerHTML = lines.join("\n");

  frames++;
  const now = performance.now();
  if (now - fpsT0 >= 500) {
    lastFps = (frames * 1000) / (now - fpsT0);
    frames = 0;
    fpsT0 = now;
  }
}

function rebuild() {
  nx = parseInt(resSel.value, 10);
  const scenario = scenarioSel.value;
  const dt = scenario === "blowup" ? 2e-3 : 5e-3;
  lx = scenario === "akhmediev" ? Math.PI / Math.sin(1.0) : 80.0;
  sim = new NlsSim(nx, lx, dt, scenario, 42);
  yoff = 0;
  etaEl.className = "";
  etaEl.textContent = "—";
}

function frame() {
  if (!paused && sim) {
    for (let s = 0; s < subSteps; s++) sim.step();
    draw();
  }
  requestAnimationFrame(frame);
}

async function main() {
  await init();
  resize();
  window.addEventListener("resize", resize);
  rebuild();

  scenarioSel.addEventListener("change", rebuild);
  resSel.addEventListener("change", rebuild);
  document.getElementById("reset")!.addEventListener("click", rebuild);
  document.getElementById("pause")!.addEventListener("click", (e) => {
    paused = !paused;
    (e.target as HTMLButtonElement).textContent = paused ? "Resume" : "Pause";
  });
  speedEl.addEventListener("input", () => {
    subSteps = parseInt(speedEl.value, 10);
    speedVal.textContent = `${subSteps}×`;
  });

  requestAnimationFrame(frame);
}

main().catch((err) => {
  statsEl.textContent = `init failed:\n${err}`;
  console.error(err);
});
