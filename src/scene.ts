// setup scene

import { Acceleration } from "./accelerometer";
import { FlipFluid } from "./flipFluid";
export interface SceneOptions {
  simWidth: number,
  simHeight: number,
  density?: number,
}

export class Scene {
  public gravity: Acceleration;
  public dt: number;
  public flipRatio: number;
  public numPressureIters: number;
  public numParticleIters: number;
  public overRelaxation: number;
  public compensateDrift: boolean;
  public separateParticles: boolean;
  public fluid: FlipFluid;
  public simWidth: number;
  public simHeight: number;
  public spacing: number;
  public particleRadius: number;
  public dynNumX: number;
  public dynNumY: number;
  public density: number;

  constructor() {
    this.gravity = {
      x: 0,
      y: -9.80,
    };
    this.dt = 1.0 / 120.0;
    this.flipRatio = 0.9;
    this.numPressureIters = 100;
    this.numParticleIters = 2;
    this.overRelaxation = 1.9;
    this.compensateDrift = true;
    this.separateParticles = true;
    this.fluid = new FlipFluid({
      density: 0,
      simWidth: 0,
      simHeight: 0,
      spacing: 0,
      particleRadius: 0,
      dynNumX: 100,
      dynNumY: 100,
    });
    this.simWidth = 0;
    this.simHeight = 0;
    this.spacing = 0;
    this.particleRadius = 0;
    this.dynNumX = 100;
    this.dynNumY = 100;
    this.density = 0;
  }


  setupScene(opts: SceneOptions) {
    const { simWidth, simHeight } = opts;

    this.simHeight = simHeight;
    this.simWidth = simWidth;


    this.overRelaxation = 1.9;

    this.dt = 1.0 / 60.0;
    this.numPressureIters = 50;
    this.numParticleIters = 2;

    const res = 100;

    const spacing = Math.max(simHeight, simWidth) / res;
    this.spacing = spacing;
    const density = opts.density ?? 1000.0;
    this.density = density;
    // compute number of particles
    const particleRadius = 0.3 * spacing;  // particle radius w.r.t. cell size
    this.particleRadius = particleRadius;
    // create fluid
    this.fluid = new FlipFluid({
      density, simWidth, simHeight, spacing, particleRadius,
      dynNumX: 100, dynNumY: 100
    });
  }
}
