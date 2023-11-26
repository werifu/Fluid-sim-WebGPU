async function main() {
	const app = document.querySelector<HTMLDivElement>('#app')!;
	const canvas = document.querySelector<HTMLCanvasElement>('#container')!;

	//   canvas.width = window.innerWidth - 20;
	//   canvas.height = window.innerHeight - 100;
	canvas.focus();
	const simHeight = 3.0;
	const cScale = canvas.height / simHeight;
	const simWidth = canvas.width / cScale;

	// const U_FIELD = 0;
	// const V_FIELD = 1;

	const FLUID_CELL = 0;
	const AIR_CELL = 1;
	const SOLID_CELL = 2;

	// const cnt = 0;

	const gl = canvas.getContext('webgl')!;
	// const vertexShaderSource = `
	//   attribute vec2 a_position;
	//   uniform vec2 u_resolution;
	//   uniform vec2 u_translation;
	//   void main() {
	//     vec2 position = a_position + u_translation;
	//     vec2 clipSpace = (position / u_resolution) * 2.0 - 1.0;
	//     gl_Position = vec4(clipSpace * vec2(1, -1), 0, 1);
	//   }
	// `;

	// // control the precision of float
	// // lowp: 8 bits float
	// // mediump: 16 bits float
	// // highp: 32 bits float
	// const fragmentShaderSource = `
	//   precision mediump float;
	//   uniform vec4 u_color;
	//   void main() {
	//     gl_FragColor = u_color;
	//   }
	// `;

	// shade create helper for vertex shader and fragement shader
	function loadShader(gl: WebGLRenderingContext, type: number, source: WebGLShader & string) {
		// "!" makes sure 
		const shader = gl.createShader(type)!;
		gl.shaderSource(shader, source);
		gl.compileShader(shader);
		// check the compile status of the shader
		// if not true
		if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
			// log out the error
			console.error('An error occurred compiling the shaders:', gl.getShaderInfoLog(shader));
			// remove the shader
			gl.deleteShader(shader);
			return null;
		}
		return shader;
	}

	function clamp(x: number, min: number, max: number) {
		if (x < min)
			return min;
		else if (x > max)
			return max;
		else
			return x;
	}

	class FlipFluid {
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

		constructor(density: number, width: number, height: number, spacing: number, particleRadius: number, maxParticles: number) {
			// fluid

			this.density = density;

			// distance between central point of one cell and neighbour
			// 600 / 6 + 1 = 101
			this.fNumX = Math.floor(width / spacing) + 1;
			// 300 / 6 + 1 = 51
			this.fNumY = Math.floor(height / spacing) + 1;
			// height(width) of each cell
			this.h = Math.max(width / this.fNumX, height / this.fNumY);
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

			// particles

			// maximum number of particles
			this.maxParticles = maxParticles;

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
			this.pNumX = Math.floor(width * this.pInvSpacing) + 1;
			this.pNumY = Math.floor(height * this.pInvSpacing) + 1;
			this.pNumCells = this.pNumX * this.pNumY;

			this.numCellParticles = new Int32Array(this.pNumCells);
			this.firstCellParticle = new Int32Array(this.pNumCells + 1);
			this.cellParticleIds = new Int32Array(maxParticles);

			// total number of particles
			this.numParticles = 0;
		}

		// update velocity and location of each particle
		public integrateParticles(dt: number, gravityAccer: number) {
			for (var i = 0; i < this.numParticles; i++) {
				// velocity along y-axis: v = v + g * dt
				this.particleVel[2 * i + 1] += dt * gravityAccer;
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

		public handleParticleCollisions(obstacleX: number, obstacleY: number, obstacleRadius: number) {
			var h = 1.0 / this.fInvSpacing;
			var r = this.particleRadius;
			// var or = obstacleRadius;
			// var or2 = or * or;
			var minDist = obstacleRadius + r;
			var minDist2 = minDist * minDist;

			var minX = h + r;
			var maxX = (this.fNumX - 1) * h - r;
			var minY = h + r;
			var maxY = (this.fNumY - 1) * h - r;


			for (var i = 0; i < this.numParticles; i++) {
				var x = this.particlePos[2 * i];
				var y = this.particlePos[2 * i + 1];

				var dx = x - obstacleX;
				var dy = y - obstacleY;
				var d2 = dx * dx + dy * dy;

				// obstacle collision

				if (d2 < minDist2) {
					this.particleVel[2 * i] = scene.obstacleVelX;
					this.particleVel[2 * i + 1] = scene.obstacleVelY;
				}

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

		updateParticleDensity(f: FlipFluid) {
			var n = this.fNumY;
			var h = this.h;
			var h1 = this.fInvSpacing;
			var h2 = 0.5 * h;

			var d = f.particleDensity;

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

		simulate(dt: number, gravity: number, flipRatio: number, numPressureIters: number, numParticleIters: number,
			overRelaxation: number, compensateDrift: boolean, separateParticles: boolean,
			obstacleX: number, abstacleY: number, obstacleRadius: number, f: FlipFluid) {
			var numSubSteps = 1;
			var sdt = dt / numSubSteps;

			for (var step = 0; step < numSubSteps; step++) {
				this.integrateParticles(sdt, gravity);
				if (separateParticles)
					this.pushParticlesApart(numParticleIters);
				this.handleParticleCollisions(obstacleX, abstacleY, obstacleRadius)
				this.transferVelocities(true, flipRatio);
				this.updateParticleDensity(f);
				this.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
				this.transferVelocities(false, flipRatio);
			}

			this.updateParticleColors();
			this.updateCellColors();

		}
	}


	const scene =
	{
		gravity: -9.80,
		dt: 1.0 / 120.0,
		flipRatio: 0.9,
		numPressureIters: 100,
		numParticleIters: 2,
		frameNr: 0,
		overRelaxation: 1.9,
		compensateDrift: true,
		separateParticles: true,
		obstacleX: 0.0,
		obstacleY: 0.0,
		obstacleRadius: canvas.height / 900,
		paused: true,
		showObstacle: true,
		obstacleVelX: 0.0,
		obstacleVelY: 0.0,
		showParticles: true,
		showGrid: false,
		fluid: new FlipFluid(0, 0, 0, 0, 0, 0)
	};
	// var f = scene.fluid = new FlipFluid(0, 0, 0, 0, 0, 0);

	function setupScene() {
		scene.obstacleRadius = 0.15;
		scene.overRelaxation = 1.9;

		scene.dt = 1.0 / 60.0;
		scene.numPressureIters = 50;
		scene.numParticleIters = 2;

		var res = 100;

		/* adapt to mobile device */
		var tankHeight = 1.0 * simHeight;
		var tankWidth = 1.0 * simWidth;
		var h = tankHeight / res;
		var density = 1000.0;

		var relWaterHeight = 0.8
		var relWaterWidth = 0.6

		// compute number of particles

		var r = 0.3 * h;  // particle radius w.r.t. cell size
		var dx = 2.0 * r;
		var dy = Math.sqrt(3.0) / 2.0 * dx;

		var numX = Math.floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
		var numY = Math.floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
		const maxParticles = numX * numY;

		// create fluid

		var f = scene.fluid = new FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

		// create particles

		f.numParticles = numX * numY;
		var p = 0;
		for (var i = 0; i < numX; i++) {
			for (var j = 0; j < numY; j++) {
				f.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
				f.particlePos[p++] = h + r + dy * j
			}
		}

		// setup grid cells for tank

		var n = f.fNumY;

		for (var i = 0; i < f.fNumX; i++) {
			for (var j = 0; j < f.fNumY; j++) {
				var s = 1.0;  // fluid
				if (i == 0 || i == f.fNumX - 1 || j == 0)
					s = 0.0;  // solid
				f.s[i * n + j] = s
			}
		}

		setObstacle(3.0, 2.6, true);
		return f;
	}

	// draw -------------------------------------------------------

	const pointVertexShader = `
    attribute vec2 attrPosition;
    attribute vec3 attrColor;
    uniform vec2 domainSize;
    uniform float pointSize;
    uniform float drawDisk;

    varying vec3 fragColor;
    varying float fragDrawDisk;

    void main() {
    vec4 screenTransform = 
    vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);
    gl_Position =
      vec4(attrPosition * screenTransform.xy + screenTransform.zw, 0.0, 1.0);

    gl_PointSize = pointSize;
    fragColor = attrColor;
    fragDrawDisk = drawDisk;
    }
  `;

	const pointFragmentShader = `
  precision mediump float;
  varying vec3 fragColor;
  varying float fragDrawDisk;

  void main() {
    if (fragDrawDisk == 1.0) {
      float rx = 0.5 - gl_PointCoord.x;
      float ry = 0.5 - gl_PointCoord.y;
      float r2 = rx * rx + ry * ry;
      if (r2 > 0.25)
        discard;
    }
    gl_FragColor = vec4(fragColor, 1.0);
  }
`;

	const meshVertexShader = `
  attribute vec2 attrPosition;
  uniform vec2 domainSize;
  uniform vec3 color;
  uniform vec2 translation;
  uniform float scale;

  varying vec3 fragColor;

  void main() {
    vec2 v = translation + attrPosition * scale;
  vec4 screenTransform = 
    vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);
  gl_Position =
    vec4(v * screenTransform.xy + screenTransform.zw, 0.0, 1.0);

  fragColor = color;
  }
`;

	const meshFragmentShader = `
  precision mediump float;
  varying vec3 fragColor;

  void main() {
    gl_FragColor = vec4(fragColor, 1.0);
  }
`;

	function createShader(gl: WebGLRenderingContext, vsSource: WebGLShader & string, fsSource: WebGLShader & string) {
		const vsShader = gl.createShader(gl.VERTEX_SHADER)!;
		gl.shaderSource(vsShader, vsSource);
		gl.compileShader(vsShader);
		if (!gl.getShaderParameter(vsShader, gl.COMPILE_STATUS))
			console.log("vertex shader compile error: " + gl.getShaderInfoLog(vsShader));

		const fsShader = gl.createShader(gl.FRAGMENT_SHADER)!;
		gl.shaderSource(fsShader, fsSource);
		gl.compileShader(fsShader);
		if (!gl.getShaderParameter(fsShader, gl.COMPILE_STATUS))
			console.log("fragment shader compile error: " + gl.getShaderInfoLog(fsShader));

		var shader = gl.createProgram()!;
		gl.attachShader(shader, vsShader);
		gl.attachShader(shader, fsShader);
		gl.linkProgram(shader);

		return shader;
	}

	let pointShader: any = null;
	var meshShader: any = null;

	var pointVertexBuffer: any = null;
	var pointColorBuffer: any = null;

	var gridVertBuffer: any = null;
	var gridColorBuffer: any = null;

	var diskVertBuffer: any = null;
	var diskIdBuffer: any = null;

	function draw() {
		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT);

		gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

		// prepare shaders

		if (pointShader == null)
			pointShader = createShader(gl, pointVertexShader, pointFragmentShader);
		if (meshShader == null)
			meshShader = createShader(gl, meshVertexShader, meshFragmentShader);

		// grid

		if (gridVertBuffer == null) {

			var f = scene.fluid;
			gridVertBuffer = gl.createBuffer();
			var cellCenters = new Float32Array(2 * f.fNumCells);
			var p = 0;

			for (var i = 0; i < f.fNumX; i++) {
				for (var j = 0; j < f.fNumY; j++) {
					cellCenters[p++] = (i + 0.5) * f.h;
					cellCenters[p++] = (j + 0.5) * f.h;
				}
			}
			gl.bindBuffer(gl.ARRAY_BUFFER, gridVertBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, cellCenters, gl.DYNAMIC_DRAW);
			gl.bindBuffer(gl.ARRAY_BUFFER, null);
		}

		if (gridColorBuffer == null)
			gridColorBuffer = gl.createBuffer();

		if (scene.showGrid) {

			var pointSize = 0.9 * scene.fluid.h / simWidth * canvas.width;

			gl.useProgram(pointShader);
			gl.uniform2f(gl.getUniformLocation(pointShader, 'domainSize'), simWidth, simHeight);
			gl.uniform1f(gl.getUniformLocation(pointShader, 'pointSize'), pointSize);
			gl.uniform1f(gl.getUniformLocation(pointShader, 'drawDisk'), 0.0);

			gl.bindBuffer(gl.ARRAY_BUFFER, gridVertBuffer);
			var posLoc = gl.getAttribLocation(pointShader, 'attrPosition');
			gl.enableVertexAttribArray(posLoc);
			gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

			gl.bindBuffer(gl.ARRAY_BUFFER, gridColorBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.cellColor, gl.DYNAMIC_DRAW);

			var colorLoc = gl.getAttribLocation(pointShader, 'attrColor');
			gl.enableVertexAttribArray(colorLoc);
			gl.vertexAttribPointer(colorLoc, 3, gl.FLOAT, false, 0, 0);

			gl.drawArrays(gl.POINTS, 0, scene.fluid.fNumCells);

			gl.disableVertexAttribArray(posLoc);
			gl.disableVertexAttribArray(colorLoc);

			gl.bindBuffer(gl.ARRAY_BUFFER, null);
		}

		// water

		if (scene.showParticles) {
			gl.clear(gl.DEPTH_BUFFER_BIT);

			var pointSize = 2.0 * scene.fluid.particleRadius / simWidth * canvas.width;

			gl.useProgram(pointShader);
			gl.uniform2f(gl.getUniformLocation(pointShader, 'domainSize'), simWidth, simHeight);
			gl.uniform1f(gl.getUniformLocation(pointShader, 'pointSize'), pointSize);
			gl.uniform1f(gl.getUniformLocation(pointShader, 'drawDisk'), 1.0);

			if (pointVertexBuffer == null)
				pointVertexBuffer = gl.createBuffer();
			if (pointColorBuffer == null)
				pointColorBuffer = gl.createBuffer();

			gl.bindBuffer(gl.ARRAY_BUFFER, pointVertexBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.particlePos, gl.DYNAMIC_DRAW);

			var posLoc = gl.getAttribLocation(pointShader, 'attrPosition');
			gl.enableVertexAttribArray(posLoc);
			gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

			gl.bindBuffer(gl.ARRAY_BUFFER, pointColorBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.particleColor, gl.DYNAMIC_DRAW);

			var colorLoc = gl.getAttribLocation(pointShader, 'attrColor');
			gl.enableVertexAttribArray(colorLoc);
			gl.vertexAttribPointer(colorLoc, 3, gl.FLOAT, false, 0, 0);

			gl.drawArrays(gl.POINTS, 0, scene.fluid.numParticles);

			gl.disableVertexAttribArray(posLoc);
			gl.disableVertexAttribArray(colorLoc);

			gl.bindBuffer(gl.ARRAY_BUFFER, null);
		}

		// disk

		// prepare disk mesh

		var numSegs = 50;

		if (diskVertBuffer == null) {

			diskVertBuffer = gl.createBuffer();
			var dphi = 2.0 * Math.PI / numSegs;
			var diskVerts = new Float32Array(2 * numSegs + 2);
			var p = 0;
			diskVerts[p++] = 0.0;
			diskVerts[p++] = 0.0;
			for (var i = 0; i < numSegs; i++) {
				diskVerts[p++] = Math.cos(i * dphi);
				diskVerts[p++] = Math.sin(i * dphi);
			}
			gl.bindBuffer(gl.ARRAY_BUFFER, diskVertBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, diskVerts, gl.DYNAMIC_DRAW);
			gl.bindBuffer(gl.ARRAY_BUFFER, null);

			diskIdBuffer = gl.createBuffer();
			var diskIds = new Uint16Array(3 * numSegs);
			p = 0;
			for (var i = 0; i < numSegs; i++) {
				diskIds[p++] = 0;
				diskIds[p++] = 1 + i;
				diskIds[p++] = 1 + (i + 1) % numSegs;
			}

			gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, diskIdBuffer);
			gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, diskIds, gl.DYNAMIC_DRAW);
			gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
		}

		gl.clear(gl.DEPTH_BUFFER_BIT);

		var diskColor = [1.0, 0.0, 0.0];

		gl.useProgram(meshShader);
		gl.uniform2f(gl.getUniformLocation(meshShader, 'domainSize'), simWidth, simHeight);
		gl.uniform3f(gl.getUniformLocation(meshShader, 'color'), diskColor[0], diskColor[1], diskColor[2]);
		gl.uniform2f(gl.getUniformLocation(meshShader, 'translation'), scene.obstacleX, scene.obstacleY);
		gl.uniform1f(gl.getUniformLocation(meshShader, 'scale'), scene.obstacleRadius + scene.fluid.particleRadius);

		posLoc = gl.getAttribLocation(meshShader, 'attrPosition');
		gl.enableVertexAttribArray(posLoc);
		gl.bindBuffer(gl.ARRAY_BUFFER, diskVertBuffer);
		gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, diskIdBuffer);
		gl.drawElements(gl.TRIANGLES, 3 * numSegs, gl.UNSIGNED_SHORT, 0);

		gl.disableVertexAttribArray(posLoc);

	}

	function setObstacle(x: number, y: number, reset: boolean) {

		var vx = 0.0;
		var vy = 0.0;

		if (!reset) {
			vx = (x - scene.obstacleX) / scene.dt;
			vy = (y - scene.obstacleY) / scene.dt;
		}

		scene.obstacleX = x;
		scene.obstacleY = y;
		var r = scene.obstacleRadius;
		var f = scene.fluid;
		var n = f.fNumY;
		// var cd = Math.sqrt(2) * f.h;

		for (var i = 1; i < f.fNumX - 2; i++) {
			for (var j = 1; j < f.fNumY - 2; j++) {

				f.s[i * n + j] = 1.0;

				var dx = (i + 0.5) * f.h - x;
				var dy = (j + 0.5) * f.h - y;

				if (dx * dx + dy * dy < r * r) {
					f.s[i * n + j] = 0.0;
					f.u[i * n + j] = vx;
					f.u[(i + 1) * n + j] = vx;
					f.v[i * n + j] = vy;
					f.v[i * n + j + 1] = vy;
				}
			}
		}

		scene.showObstacle = true;
		scene.obstacleVelX = vx;
		scene.obstacleVelY = vy;
	}

	// interaction -------------------------------------------------------

	var mouseDown = false;

	function startDrag(x: number, y: number) {
		// console.log(" start drag here");
		let bounds = canvas.getBoundingClientRect();

		let mx = x - bounds.left - canvas.clientLeft;
		let my = y - bounds.top - canvas.clientTop;
		mouseDown = true;

		x = mx / cScale;
		y = (canvas.height - my) / cScale;

		setObstacle(x, y, true);
		scene.paused = false;
	}

	function drag(x: number, y: number) {
		if (mouseDown) {
			// console.log(" drag here");
			let bounds = canvas.getBoundingClientRect();
			let mx = x - bounds.left - canvas.clientLeft;
			let my = y - bounds.top - canvas.clientTop;
			x = mx / cScale;
			y = (canvas.height - my) / cScale;
			setObstacle(x, y, false);
		}
	}

	function endDrag() {
		// console.log("end drag here");
		mouseDown = false;
		scene.obstacleVelX = 0.0;
		scene.obstacleVelY = 0.0;
	}

	canvas.addEventListener('mousedown', event => {
		// console.log("mouse down here");
		startDrag(event.x, event.y);
	});

	canvas.addEventListener('mouseup', event => {
		// console.log("mouse up here");
		endDrag();
	});

	canvas.addEventListener('mousemove', event => {
		// console.log("mousemove here");
		drag(event.x, event.y);
		// console.log("corrdinate: ", event.x, event.y);
	});

	canvas.addEventListener('touchstart', event => {
		// console.log("touchstart here");
		startDrag(event.touches[0].clientX, event.touches[0].clientY)
	});

	canvas.addEventListener('touchend', event => {
		// console.log("touchend here");
		endDrag()
	});

	canvas.addEventListener('touchmove', event => {
		// console.log("touchmove here");
		event.preventDefault();
		event.stopImmediatePropagation();
		drag(event.touches[0].clientX, event.touches[0].clientY)
	}, { passive: false });

	document.addEventListener('keydown', event => {
		switch (event.key) {
			case 'p': scene.paused = !scene.paused; break;
			case 'm': scene.paused = false; simulate(); scene.paused = true; break;
		}
	});

	let acl = {
		x: 0,
		y: scene.gravity,
	};

	window.addEventListener('devicemotion', (event: DeviceMotionEvent) => {
		acl.x = -event.accelerationIncludingGravity?.x!;
		acl.y = event.accelerationIncludingGravity?.y!;
		app.innerText = `aclX: ${acl.x}\naclY:${acl.y}`;
	})

	let totalTime = 0.0;
	var current_time = 0.0;
	const variableChangeEvent = new Event('variableChange');

	function changeVariableValue(newValue: number) {
		if (totalTime != newValue) {
			totalTime = newValue;
			canvas.dispatchEvent(variableChangeEvent);
		}
	}

	canvas.addEventListener('variableChange', function (event) {
		var vx = 0.0;
		var vy = 0.0;

		vx += (acl.x * 200) * scene.dt;
		vy += (acl.y * 200) * scene.dt;

		checkBoundary(scene.obstacleX + vx * scene.dt, scene.obstacleY - vy * scene.dt);

		var r = scene.obstacleRadius;
		var f = scene.fluid;
		var n = f.fNumY;
		// var cd = Math.sqrt(2) * f.h;

		for (var i = 1; i < f.fNumX - 2; i++) {
			for (var j = 1; j < f.fNumY - 2; j++) {

				f.s[i * n + j] = 1.0;

				var dx = (i + 0.5) * f.h - scene.obstacleX;
				var dy = (j + 0.5) * f.h - scene.obstacleY;

				if (dx * dx + dy * dy < r * r) {
					f.s[i * n + j] = 0.0;
					f.u[i * n + j] = vx;
					f.u[(i + 1) * n + j] = vx;
					f.v[i * n + j] = vy;
					f.v[i * n + j + 1] = vy;
				}
			}
		}

		scene.showObstacle = true;
		scene.obstacleVelX = vx;
		scene.obstacleVelY = vy;
	});

	function addTime() {
		if (scene.paused)
			scene.paused = false;
		current_time = totalTime;
		current_time += scene.dt;
		changeVariableValue(current_time);
	}

	setInterval(addTime, scene.dt);

	function checkBoundary(x: number, y: number) {
		// 0 3
		// 4 3
		// 0 0 
		// 4 0
		// let bounds = canvas.getBoundingClientRect();
		// let top = bounds.top;
		// let left = bounds.left;
		// let right = bounds.right;

		if ((x < 0.2) && (y >= 0.2 && y <= 2.8)) {
			scene.obstacleX = 0.2;
			scene.obstacleY = y;
		}

		if ((x > 3.8) && (y >= 0.2 && y <= 2.8)) {
			scene.obstacleX = 3.8;
			scene.obstacleY = y;
		}

		if ((y < 0.2) && (x >= 0.2 || x <= 3.8)) {
			scene.obstacleX = x;
			scene.obstacleY = 0.2;
		}

		if ((y > 2.8) && (x >= 0.2 || x <= 3.8)) {
			scene.obstacleX = x;
			scene.obstacleY = 2.8;
		}

		if (x < 0.2 && y < 0.2) {
			scene.obstacleX = 0.2;
			scene.obstacleY = 0.2;
		}
		else if (x < 0.2 && y > 2.8) {
			scene.obstacleX = 0.2;
			scene.obstacleY = 2.8;
		}
		else if (x > 3.8 && y < 0.2) {
			scene.obstacleX = 3.8;
			scene.obstacleY = 0.2;
		}
		else if (x > 3.8 && y > 2.8) {
			scene.obstacleX = 3.8;
			scene.obstacleY = 2.8;
		}

		if (x >= 0 && x <= 4 && y >= 0 && y <= 3) {
			scene.obstacleX = x;
			scene.obstacleY = y;
		}

	}

	// function toggleStart()
	// {
	//   const button = document.getElementById('startButton');
	//   if (scene.paused)
	//     button!.innerHTML = "Stop";
	//   else
	//     button!.innerHTML = "Start";
	//   scene.paused = !scene.paused;
	// }

	function simulate() {
		if (!scene.paused)
			scene.fluid.simulate(
				scene.dt, scene.gravity, scene.flipRatio, scene.numPressureIters, scene.numParticleIters,
				scene.overRelaxation, scene.compensateDrift, scene.separateParticles,
				scene.obstacleX, scene.obstacleY, scene.obstacleRadius, f);
		scene.frameNr++;
	}

	function update() {
		simulate();
		draw();
		requestAnimationFrame(update);
	}

	const f = setupScene();
	update();
}

main();