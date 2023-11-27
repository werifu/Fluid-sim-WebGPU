import { Scene } from "./scene";

interface FluidOptions {
  density: number,
  simWidth: number,
  simHeight: number,
  spacing: number,
  particleRadius: number,
  dynNumX: number,
  dynNumY: number,
}

const FLUID_CELL = 0;
const AIR_CELL = 1;
const SOLID_CELL = 2;

function clamp(x: number, min: number, max: number) {
  if (x < min)
    return min;
  else if (x > max)
    return max;
  else
    return x;
}

export class FlipFluid {
  public density: number;
  public fNumX!: number;
  public fNumY!: number;
  // grid unit
  public h!: number;
  public fInvSpacing: number;
  public fNumCells!: number;
  // velocity field: (u, v)
  public u!: Float32Array;
  public v!: Float32Array;
  public du!: Float32Array;
  public dv!: Float32Array;
  public prevU!: Float32Array;
  public prevV!: Float32Array;
  public p!: Float32Array;
  public s!: Float32Array;
  // cell type: water or air
  public cellType!: Int32Array;
  public cellColor!: Float32Array;

  // maximum number of particles
  public maxParticles!: number;
  // particles position list
  public particlePos!: Float32Array;
  // particles color list
  public particleColor: Float32Array;
  // particles velocity list
  public particleVel: Float32Array;
  // particles density list
  public particleDensity: Float32Array;
  // the rest of particles density list
  public particleRestDensity: number;
  // radius of particle
  public particleRadius: number;
  public pInvSpacing: number;
  public pNumX: number;
  public pNumY: number;
  public pNumCells: number;
  public numCellParticles: Int32Array;
  public firstCellParticle: Int32Array;
  public cellParticleIds: Int32Array;
  // total number of particles
  public numParticles: number;

  constructor(opts: FluidOptions) {
    const {
      density,
      simWidth,
      simHeight,
      spacing,
      particleRadius,
      dynNumX,
      dynNumY,
    } = opts;

    // total number of particles

    const res = 100;

    const h = Math.max(simHeight, simWidth) / res;

    // compute number of particles
    // const relWaterHeight = 0.8;
    // const relWaterWidth = 0.6;
    const r = 0.3 * h;  // particle radius w.r.t. cell size
    const dx = 2.0 * r;
    const dy = Math.sqrt(3.0) / 2.0 * dx;

    // const numX = Math.min(dynNumX, Math.floor((relWaterWidth * simWidth - 2.0 * h - 2.0 * r) / dx));
    // const numY = Math.min(dynNumY, Math.floor((relWaterHeight * simHeight - 2.0 * h - 2.0 * r) / dy));
    const numX = dynNumX;
    const numY = dynNumY;
    // const maxParticles = numX * numY;

    this.numParticles = numX * numY;
    // maximum number of particles
    this.maxParticles = this.numParticles + 1;

    // fluid
    this.density = density;

    // distance between central point of one cell and neighbour
    // 600 / 6 + 1 = 101
    this.fNumX = Math.floor(simWidth / spacing) + 1;
    // 300 / 6 + 1 = 51
    this.fNumY = Math.floor(simHeight / spacing) + 1;
    // height(width) of each cell
    this.h = Math.max(simWidth / this.fNumX, simHeight / this.fNumY);
    this.fInvSpacing = 1.0 / this.h;
    this.fNumCells = this.fNumX * this.fNumY;

    this.u = new Float32Array(this.fNumCells);
    this.v = new Float32Array(this.fNumCells);
    this.du = new Float32Array(this.fNumCells);
    this.dv = new Float32Array(this.fNumCells);
    this.prevU = new Float32Array(this.fNumCells);
    this.prevV = new Float32Array(this.fNumCells);
    this.p = new Float32Array(this.fNumCells);
    this.s = new Float32Array(this.fNumCells);
    this.cellType = new Int32Array(this.fNumCells);
    this.cellColor = new Float32Array(3 * this.fNumCells);


    // positions of each particle (x, y)
    this.particlePos = new Float32Array(2 * this.maxParticles);
    // color of each particle (r, g, b)
    this.particleColor = new Float32Array(3 * this.maxParticles);
    for (var i = 0; i < this.maxParticles; i++) {
      //  (r, g, b): b = 1.0 -> dark blue
      this.particleColor[3 * i + 2] = 1.0;
    }

    // velocity of each particle (u, v)
    this.particleVel = new Float32Array(2 * this.maxParticles);
    // density of particles in each cell
    this.particleDensity = new Float32Array(this.fNumCells);
    this.particleRestDensity = 0.0;

    this.particleRadius = particleRadius;
    this.pInvSpacing = 1.0 / (2.2 * particleRadius);
    this.pNumX = Math.floor(simWidth * this.pInvSpacing) + 1;
    this.pNumY = Math.floor(simHeight * this.pInvSpacing) + 1;
    this.pNumCells = this.pNumX * this.pNumY;

    this.numCellParticles = new Int32Array(this.pNumCells);
    this.firstCellParticle = new Int32Array(this.pNumCells + 1);
    this.cellParticleIds = new Int32Array(this.maxParticles);


    let p = 0;
    for (let i = 0; i < numX; i++) {
      for (let j = 0; j < numY; j++) {
        this.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
        this.particlePos[p++] = h + r + dy * j
      }
    }

    // setup grid cells for tank
    var n = this.fNumY;

    for (let i = 0; i < this.fNumX; i++) {
      for (let j = 0; j < this.fNumY; j++) {
        let s = 1.0;  // fluid
        if (i == 0 || i == this.fNumX - 1 || j == 0)
          s = 0.0;  // solid
        this.s[i * n + j] = s
      }
    }
  }

  // update velocity and location of each particle
  public integrateParticles(dt: number, gravityAccer: { x: number, y: number }) {
    for (var i = 0; i < this.numParticles; i++) {
      // velocity along x-axis: v = v + g * dt
      this.particleVel[2 * i] += dt * gravityAccer.x;
      // velocity along y-axis: v = v + g * dt
      this.particleVel[2 * i + 1] += dt * gravityAccer.y;
      // x = x + Vx * dt
      this.particlePos[2 * i] += this.particleVel[2 * i] * dt;
      // y = y + Vy * dt
      this.particlePos[2 * i + 1] += this.particleVel[2 * i + 1] * dt;
    }
  }

  public pushParticlesApart(numIters: number) {
    var colorDiffusionCoeff = 0.001;

    // count particles per cell

    this.numCellParticles.fill(0);

    for (var i = 0; i < this.numParticles; i++) {
      var x = this.particlePos[2 * i];
      var y = this.particlePos[2 * i + 1];

      var xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      var yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      var cellNr = xi * this.pNumY + yi;
      this.numCellParticles[cellNr]++;
    }

    // partial sums

    var first = 0;

    for (var i = 0; i < this.pNumCells; i++) {
      first += this.numCellParticles[i];
      this.firstCellParticle[i] = first;
    }
    this.firstCellParticle[this.pNumCells] = first;    // guard

    // fill particles into cells

    for (var i = 0; i < this.numParticles; i++) {
      var x = this.particlePos[2 * i];
      var y = this.particlePos[2 * i + 1];

      var xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      var yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      var cellNr = xi * this.pNumY + yi;
      this.firstCellParticle[cellNr]--;
      this.cellParticleIds[this.firstCellParticle[cellNr]] = i;
    }

    // push particles apart

    var minDist = 2.0 * this.particleRadius;
    var minDist2 = minDist * minDist;

    for (var iter = 0; iter < numIters; iter++) {

      for (var i = 0; i < this.numParticles; i++) {
        var px = this.particlePos[2 * i];
        var py = this.particlePos[2 * i + 1];

        var pxi = Math.floor(px * this.pInvSpacing);
        var pyi = Math.floor(py * this.pInvSpacing);
        var x0 = Math.max(pxi - 1, 0);
        var y0 = Math.max(pyi - 1, 0);
        var x1 = Math.min(pxi + 1, this.pNumX - 1);
        var y1 = Math.min(pyi + 1, this.pNumY - 1);

        for (var xi = x0; xi <= x1; xi++) {
          for (var yi = y0; yi <= y1; yi++) {
            var cellNr = xi * this.pNumY + yi;
            var first = this.firstCellParticle[cellNr];
            var last = this.firstCellParticle[cellNr + 1];
            for (var j = first; j < last; j++) {
              var id = this.cellParticleIds[j];
              if (id == i)
                continue;
              var qx = this.particlePos[2 * id];
              var qy = this.particlePos[2 * id + 1];

              var dx = qx - px;
              var dy = qy - py;
              var d2 = dx * dx + dy * dy;
              if (d2 > minDist2 || d2 == 0.0)
                continue;
              var d = Math.sqrt(d2);
              var s = 0.5 * (minDist - d) / d;
              dx *= s;
              dy *= s;

              /* change position */
              this.particlePos[2 * i] -= dx;
              this.particlePos[2 * i + 1] -= dy;
              this.particlePos[2 * id] += dx;
              this.particlePos[2 * id + 1] += dy;

              // diffuse colors

              for (var k = 0; k < 3; k++) {
                var color0 = this.particleColor[3 * i + k];
                var color1 = this.particleColor[3 * id + k];
                var color = (color0 + color1) * 0.5;
                this.particleColor[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
                this.particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
              }
            }
          }
        }
      }
    }
  }

  public handleParticleCollisions() {
    var h = 1.0 / this.fInvSpacing;
    var r = this.particleRadius;
    // var or = obstacleRadius;
    // var or2 = or * or;


    var minX = h + r;
    var maxX = (this.fNumX - 1) * h - r;
    var minY = h + r;
    var maxY = (this.fNumY - 1) * h - r;


    for (var i = 0; i < this.numParticles; i++) {
      var x = this.particlePos[2 * i];
      var y = this.particlePos[2 * i + 1];

      // wall collisions
      if (x < minX) {
        x = minX;
        this.particleVel[2 * i] = 0.0;

      }
      if (x > maxX) {
        x = maxX;
        this.particleVel[2 * i] = 0.0;
      }
      if (y < minY) {
        y = minY;
        this.particleVel[2 * i + 1] = 0.0;
      }
      if (y > maxY) {
        y = maxY;
        this.particleVel[2 * i + 1] = 0.0;
      }
      this.particlePos[2 * i] = x;
      this.particlePos[2 * i + 1] = y;
    }
  }

  updateParticleDensity() {
    var n = this.fNumY;
    var h = this.h;
    var h1 = this.fInvSpacing;
    var h2 = 0.5 * h;

    var d = this.particleDensity;

    d.fill(0.0);

    for (var i = 0; i < this.numParticles; i++) {
      var x = this.particlePos[2 * i];
      var y = this.particlePos[2 * i + 1];

      x = clamp(x, h, (this.fNumX - 1) * h);
      y = clamp(y, h, (this.fNumY - 1) * h);

      var x0 = Math.floor((x - h2) * h1);
      var tx = ((x - h2) - x0 * h) * h1;
      var x1 = Math.min(x0 + 1, this.fNumX - 2);

      var y0 = Math.floor((y - h2) * h1);
      var ty = ((y - h2) - y0 * h) * h1;
      var y1 = Math.min(y0 + 1, this.fNumY - 2);

      var sx = 1.0 - tx;
      var sy = 1.0 - ty;

      if (x0 < this.fNumX && y0 < this.fNumY) d[x0 * n + y0] += sx * sy;
      if (x1 < this.fNumX && y0 < this.fNumY) d[x1 * n + y0] += tx * sy;
      if (x1 < this.fNumX && y1 < this.fNumY) d[x1 * n + y1] += tx * ty;
      if (x0 < this.fNumX && y1 < this.fNumY) d[x0 * n + y1] += sx * ty;
    }

    if (this.particleRestDensity == 0.0) {
      var sum = 0.0;
      var numFluidCells = 0;

      for (var i = 0; i < this.fNumCells; i++) {
        if (this.cellType[i] == FLUID_CELL) {
          sum += d[i];
          numFluidCells++;
        }
      }

      if (numFluidCells > 0)
        this.particleRestDensity = sum / numFluidCells;
    }
  }

  transferVelocities(toGrid: boolean, flipRatio: number) {
    var n = this.fNumY;
    var h = this.h;
    var h1 = this.fInvSpacing;
    var h2 = 0.5 * h;

    if (toGrid) {

      this.prevU.set(this.u);
      this.prevV.set(this.v);

      this.du.fill(0.0);
      this.dv.fill(0.0);
      this.u.fill(0.0);
      this.v.fill(0.0);

      for (var i = 0; i < this.fNumCells; i++)
        this.cellType[i] = this.s[i] == 0.0 ? SOLID_CELL : AIR_CELL;

      for (var i = 0; i < this.numParticles; i++) {
        var x = this.particlePos[2 * i];
        var y = this.particlePos[2 * i + 1];
        var xi = clamp(Math.floor(x * h1), 0, this.fNumX - 1);
        var yi = clamp(Math.floor(y * h1), 0, this.fNumY - 1);
        var cellNr = xi * n + yi;
        if (this.cellType[cellNr] == AIR_CELL)
          this.cellType[cellNr] = FLUID_CELL;
      }
    }

    for (var component = 0; component < 2; component++) {

      var dx = component == 0 ? 0.0 : h2;
      var dy = component == 0 ? h2 : 0.0;

      var f = component == 0 ? this.u : this.v;
      var prevF = (component == 0) ? this.prevU : this.prevV;
      var d = component == 0 ? this.du : this.dv;

      for (var i = 0; i < this.numParticles; i++) {
        var x = this.particlePos[2 * i];
        var y = this.particlePos[2 * i + 1];

        x = clamp(x, h, (this.fNumX - 1) * h);
        y = clamp(y, h, (this.fNumY - 1) * h);

        var x0 = Math.min(Math.floor((x - dx) * h1), this.fNumX - 2);
        var tx = ((x - dx) - x0 * h) * h1;
        var x1 = Math.min(x0 + 1, this.fNumX - 2);

        var y0 = Math.min(Math.floor((y - dy) * h1), this.fNumY - 2);
        var ty = ((y - dy) - y0 * h) * h1;
        var y1 = Math.min(y0 + 1, this.fNumY - 2);

        var sx = 1.0 - tx;
        var sy = 1.0 - ty;

        var d0 = sx * sy;
        var d1 = tx * sy;
        var d2 = tx * ty;
        var d3 = sx * ty;

        var nr0 = x0 * n + y0;
        var nr1 = x1 * n + y0;
        var nr2 = x1 * n + y1;
        var nr3 = x0 * n + y1;

        if (toGrid) {
          var pv = this.particleVel[2 * i + component];
          f[nr0] += pv * d0; d[nr0] += d0;
          f[nr1] += pv * d1; d[nr1] += d1;
          f[nr2] += pv * d2; d[nr2] += d2;
          f[nr3] += pv * d3; d[nr3] += d3;
        }
        else {
          var offset = component == 0 ? n : 1;
          var valid0 = this.cellType[nr0] != AIR_CELL || this.cellType[nr0 - offset] != AIR_CELL ? 1.0 : 0.0;
          var valid1 = this.cellType[nr1] != AIR_CELL || this.cellType[nr1 - offset] != AIR_CELL ? 1.0 : 0.0;
          var valid2 = this.cellType[nr2] != AIR_CELL || this.cellType[nr2 - offset] != AIR_CELL ? 1.0 : 0.0;
          var valid3 = this.cellType[nr3] != AIR_CELL || this.cellType[nr3 - offset] != AIR_CELL ? 1.0 : 0.0;

          var v = this.particleVel[2 * i + component];
          const d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

          if (d > 0.0) {

            var picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
            var corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1])
              + valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
            var flipV = v + corr;

            this.particleVel[2 * i + component] = (1.0 - flipRatio) * picV + flipRatio * flipV;
          }
        }
      }

      if (toGrid) {
        for (var i = 0; i < f.length; i++) {
          if (d[i] > 0.0)
            f[i] /= d[i];
        }

        // restore solid cells

        for (var i = 0; i < this.fNumX; i++) {
          for (var j = 0; j < this.fNumY; j++) {
            var solid = this.cellType[i * n + j] == SOLID_CELL;
            if (solid || (i > 0 && this.cellType[(i - 1) * n + j] == SOLID_CELL))
              this.u[i * n + j] = this.prevU[i * n + j];
            if (solid || (j > 0 && this.cellType[i * n + j - 1] == SOLID_CELL))
              this.v[i * n + j] = this.prevV[i * n + j];
          }
        }
      }
    }
  }

  solveIncompressibility(numIters: number, dt: number, overRelaxation: number, compensateDrift = true) {

    this.p.fill(0.0);
    this.prevU.set(this.u);
    this.prevV.set(this.v);

    var n = this.fNumY;
    var cp = this.density * this.h / dt;

    // for (var i = 0; i < this.fNumCells; i++) {
    //   var u = this.u[i];
    //   var v = this.v[i];
    // }

    for (var iter = 0; iter < numIters; iter++) {

      for (var i = 1; i < this.fNumX - 1; i++) {
        for (var j = 1; j < this.fNumY - 1; j++) {

          if (this.cellType[i * n + j] != FLUID_CELL)
            continue;

          var center = i * n + j;
          var left = (i - 1) * n + j;
          var right = (i + 1) * n + j;
          var bottom = i * n + j - 1;
          var top = i * n + j + 1;

          var s = this.s[center];
          var sx0 = this.s[left];
          var sx1 = this.s[right];
          var sy0 = this.s[bottom];
          var sy1 = this.s[top];
          var s = sx0 + sx1 + sy0 + sy1;
          if (s == 0.0)
            continue;

          var div = this.u[right] - this.u[center] +
            this.v[top] - this.v[center];

          if (this.particleRestDensity > 0.0 && compensateDrift) {
            var k = 1.0;
            var compression = this.particleDensity[i * n + j] - this.particleRestDensity;
            if (compression > 0.0)
              div = div - k * compression;
          }

          var p = -div / s;
          p *= overRelaxation;
          this.p[center] += cp * p;

          this.u[center] -= sx0 * p;
          this.u[right] += sx1 * p;
          this.v[center] -= sy0 * p;
          this.v[top] += sy1 * p;
        }
      }
    }
  }

  updateParticleColors() {
    var h1 = this.fInvSpacing;

    for (var i = 0; i < this.numParticles; i++) {

      var s = 0.01;

      this.particleColor[3 * i] = clamp(this.particleColor[3 * i] - s, 0.0, 1.0);
      this.particleColor[3 * i + 1] = clamp(this.particleColor[3 * i + 1] - s, 0.0, 1.0);
      this.particleColor[3 * i + 2] = clamp(this.particleColor[3 * i + 2] + s, 0.0, 1.0);

      var x = this.particlePos[2 * i];
      var y = this.particlePos[2 * i + 1];
      var xi = clamp(Math.floor(x * h1), 1, this.fNumX - 1);
      var yi = clamp(Math.floor(y * h1), 1, this.fNumY - 1);
      var cellNr = xi * this.fNumY + yi;

      var d0 = this.particleRestDensity;

      if (d0 > 0.0) {
        var relDensity = this.particleDensity[cellNr] / d0;
        if (relDensity < 0.7) {
          var s = 0.8;
          this.particleColor[3 * i] = s;
          this.particleColor[3 * i + 1] = s;
          this.particleColor[3 * i + 2] = 1.0;
        }
      }
    }
  }

  setSciColor(cellNr: any, val: number, minVal: number, maxVal: number) {
    val = Math.min(Math.max(val, minVal), maxVal - 0.0001);
    var d = maxVal - minVal;
    val = d == 0.0 ? 0.5 : (val - minVal) / d;
    var m = 0.25;
    var num = Math.floor(val / m);
    var s = (val - num * m) / m;
    var r = 0.0, g = 0.0, b = 0.0;

    switch (num) {
      case 0: r = 0.0; g = s; b = 1.0; break;
      case 1: r = 0.0; g = 1.0; b = 1.0 - s; break;
      case 2: r = s; g = 1.0; b = 0.0; break;
      case 3: r = 1.0; g = 1.0 - s; b = 0.0; break;
    }

    this.cellColor[3 * cellNr] = r;
    this.cellColor[3 * cellNr + 1] = g;
    this.cellColor[3 * cellNr + 2] = b;
  }

  updateCellColors() {
    this.cellColor.fill(0.0);

    for (var i = 0; i < this.fNumCells; i++) {

      if (this.cellType[i] == SOLID_CELL) {
        this.cellColor[3 * i] = 0.5;
        this.cellColor[3 * i + 1] = 0.5;
        this.cellColor[3 * i + 2] = 0.5;
      }
      else if (this.cellType[i] == FLUID_CELL) {
        var d = this.particleDensity[i];
        if (this.particleRestDensity > 0.0)
          d /= this.particleRestDensity;
        this.setSciColor(i, d, 0.0, 2.0);
      }
    }
  }

  simulate(scene: Scene) {
    const { dt, gravity, flipRatio, numPressureIters, numParticleIters,
      overRelaxation, compensateDrift, separateParticles } = scene;
    const numSubSteps = 1;
    let sdt = dt / numSubSteps;

    for (let step = 0; step < numSubSteps; step++) {
      this.integrateParticles(sdt, gravity);
      if (separateParticles)
        this.pushParticlesApart(numParticleIters);
      this.handleParticleCollisions();
      this.transferVelocities(true, flipRatio);
      this.updateParticleDensity();
      this.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
      this.transferVelocities(false, flipRatio);
    }

    this.updateParticleColors();
    this.updateCellColors();

  }
}