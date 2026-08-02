import init, { sea_scale, production_forecast } from "../pkg/rogue_wasm.js";

const hsEl = document.getElementById("hs") as HTMLInputElement;
const tpEl = document.getElementById("tp") as HTMLInputElement;
const gaEl = document.getElementById("gamma") as HTMLInputElement;
const runBtn = document.getElementById("run") as HTMLButtonElement;
const cardsEl = document.getElementById("cards")!;
const eventsEl = document.getElementById("events")!;
const jsonEl = document.getElementById("json")!;
const statusEl = document.getElementById("status")!;

function num(v: string, d: number) {
  const n = parseFloat(v);
  return Number.isFinite(n) ? n : d;
}

function card(label: string, value: string, kls = "") {
  return `<div class="card ${kls}"><div class="k">${label}</div><div class="v">${value}</div></div>`;
}

function riskClass(cf: number, note: string) {
  if (note && note.includes("SEVERE")) return "sev";
  if (cf >= 2.6) return "hi";
  if (cf >= 2.2) return "elev";
  return "low";
}

function render(f: any, sc: any) {
  const cf = f.max_crest_factor;
  const rc = riskClass(cf, f.note);
  cardsEl.innerHTML =
    card("Crest factor", cf.toFixed(2) + " σ", rc) +
    card("Max amplitude", f.max_amplitude_m.toFixed(2) + " m") +
    card("Wavelength λ₀", sc.wavelength_m.toFixed(1) + " m") +
    card("K₀", sc.k0.toFixed(4) + " /m") +
    card("Dispersion β", sc.beta.toExponential(2) + " m²/s") +
    card("Nonlinearity γ", sc.gamma_c.toExponential(2) + " /m·s") +
    card("Steepness 2Ak₀", sc.steepness.toFixed(4)) +
    card("Rogue events", String((f.events && f.events.length) || 0)) +
    card("Risk", f.note || "—", rc);

  eventsEl.innerHTML = f.events && f.events.length
    ? f.events.map((e: any) =>
        `<div class="ev">t=<b>${e.t.toFixed(1)}</b>s x=<b>${e.x_m.toFixed(1)}</b>m crest=<b>${e.crest_factor.toFixed(2)}σ</b> Hs=<b>${e.hs_m.toFixed(2)}</b>m</div>`
      ).join("")
    : `<div class="muted">No resolved rogue events this run.</div>`;

  jsonEl.innerHTML = `<pre>${JSON.stringify(f, null, 2)}</pre>`;
}

async function run() {
  const hs = num(hsEl.value, 4.0);
  const tp = num(tpEl.value, 11.0);
  const gamma = num(gaEl.value, 3.3);
  statusEl.textContent = `computing forecast for Hs=${hs} m, Tp=${tp} s …`;
  await new Promise((r) => setTimeout(r, 30));
  const sc = sea_scale(hs, tp, gamma);
  const f = production_forecast(hs, tp, gamma, 42);
  render(f, sc);
  statusEl.textContent = "done — seeded, reproducible";
}

async function main() {
  await init();
  runBtn.addEventListener("click", run);
  await run();
}

main().catch((err) => {
  statusEl.textContent = `init failed: ${err}`;
  console.error(err);
});