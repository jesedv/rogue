/* tslint:disable */
/* eslint-disable */

/**
 * 1D NLS simulation wrapped for the browser dashboard.
 */
export class NlsSim {
    free(): void;
    [Symbol.dispose](): void;
    amplitudes(): Float32Array;
    blow_up_state(): any;
    diagnostics(): any;
    /**
     * `scenario`: `akhmediev | peregrine | ocean | blowup | soliton | stokes`.
     */
    constructor(nx: number, lx: number, dt: number, scenario: string, seed: number);
    rogue_stats(): any;
    spectrum(): Float32Array;
    step(): void;
    step_count(): bigint;
    step_n(n: number): void;
    surface(): Float32Array;
    time(): number;
}

export function init(): void;

export function production_forecast(hs: number, tp: number, gamma: number, seed: number): any;

/**
 * Lighter grid for the browser demo so a click returns in well under a
 * second: fewer carriers per box and a shorter integration horizon but the
 * same physical scaling math.
 */
export function production_forecast_quick(hs: number, tp: number, gamma: number, seed: number): any;

export function sea_scale(hs: number, tp: number, gamma: number): any;

export function version(): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_nlssim_free: (a: number, b: number) => void;
    readonly nlssim_amplitudes: (a: number) => [number, number];
    readonly nlssim_blow_up_state: (a: number) => any;
    readonly nlssim_diagnostics: (a: number) => any;
    readonly nlssim_new: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number];
    readonly nlssim_rogue_stats: (a: number) => any;
    readonly nlssim_spectrum: (a: number) => [number, number];
    readonly nlssim_step: (a: number) => void;
    readonly nlssim_step_count: (a: number) => bigint;
    readonly nlssim_step_n: (a: number, b: number) => void;
    readonly nlssim_surface: (a: number) => [number, number];
    readonly nlssim_time: (a: number) => number;
    readonly production_forecast_quick: (a: number, b: number, c: number, d: number) => any;
    readonly sea_scale: (a: number, b: number, c: number) => any;
    readonly version: () => [number, number];
    readonly init: () => void;
    readonly production_forecast: (a: number, b: number, c: number, d: number) => any;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
