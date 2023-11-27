import { Acceleration } from "./accelerometer";
import { setupGUI } from "./gui";
import { Scene } from "./scene";
import { createShader, meshFragmentShader, meshVertexShader, pointFragmentShader, pointVertexShader } from "./shaders";


function permissionIOS(scene: Scene, acl: Acceleration) {
  if (typeof (DeviceMotionEvent) !== "undefined" && typeof ((DeviceMotionEvent as any).requestPermission) === "function") {
    // (optional) Do something before API request prompt.
    (DeviceMotionEvent as any).requestPermission()
      .then((response: any) => {
        // (optional) Do something after API prompt dismissed.
        if (response === "granted") {
          window.addEventListener("devicemotion", (event: DeviceMotionEvent) => {
            acl.x = event.accelerationIncludingGravity?.x!;
            acl.y = event.accelerationIncludingGravity?.y!;
            const screenAngle = window.orientation ?? 0;
            if (screenAngle == 0) {
              scene.gravity = {
                x: acl.x,
                y: acl.y
              };
            } else if (screenAngle === -90) {
              scene.gravity = {
                x: acl.y,
                y: -acl.x,
              }
            } else if (screenAngle == 90) {
              scene.gravity = {
                x: -acl.y,
                y: acl.x,
              }
            }
          })
        }
      })
      .catch(console.error)
  }
}

async function main() {
  // window.scrollTo(0, 1);
  // const app = document.querySelector<HTMLDivElement>('#app')!;
  const canvas = document.querySelector<HTMLCanvasElement>('#container')!;
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;
  canvas.focus();
  let simHeight = 4.0;
  const cScale = canvas.height / simHeight;
  let simWidth = canvas.width / cScale;
  const gl = canvas.getContext('webgl')!;

  const scene = new Scene();


  // draw -------------------------------------------------------

  let pointShader: any = null;
  var meshShader: any = null;

  var pointVertexBuffer: any = null;
  var pointColorBuffer: any = null;

  var gridVertBuffer: any = null;
  var gridColorBuffer: any = null;


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

    // water

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

  // phone orientation changed
  if (screen.orientation) {
    screen.orientation.onchange = (_: Event) => {
      canvas.width = window.innerWidth;
      canvas.height = window.innerHeight;

      if (screen.orientation.angle === 90 || screen.orientation.angle === 270) {
        simWidth = 4.0;
        const cScale = canvas.width / simWidth;
        simHeight = canvas.height / cScale;
      } else {
        simHeight = 4.0;
        const cScale = canvas.height / simHeight;
        simWidth = canvas.width / cScale;
      }
      scene.setupScene({ simWidth, simHeight, density: 1000 });
    };
  } else if (window.orientation !== undefined) {
    window.addEventListener("resize", () => {
      // Get screen size (inner/outerWidth, inner/outerHeight)
      canvas.width = window.innerWidth;
      canvas.height = window.innerHeight;

      if (window.orientation === 90 || window.orientation === -90) {
        simWidth = 4.0;
        const cScale = canvas.width / simWidth;
        simHeight = canvas.height / cScale;
      } else {
        simHeight = 4.0;
        const cScale = canvas.height / simHeight;
        simWidth = canvas.width / cScale;
      }
      scene.setupScene({ simWidth, simHeight, density: 1000 });
    }, false);
  }


  let acl: Acceleration = scene.gravity;

  // for ios devices
  if (/iPad|iPhone|Macintosh/.test(navigator.userAgent)) {
    var button = document.createElement('button');
    button.innerHTML = 'iOS Click Me!';

    // Style the button
    button.style.position = 'fixed';
    button.style.fontSize = '30px';
    button.style.top = '50%';
    button.style.left = '50%';
    button.style.transform = 'translate(-50%, -50%)';

    // Add the button to the document body
    document.body.appendChild(button);

    // Handle the button click event
    button.addEventListener('click', () => {
      permissionIOS(scene, acl);
      button.remove();
    });
  } else {
    window.addEventListener('devicemotion', (event: DeviceMotionEvent) => {
      acl.x = event.accelerationIncludingGravity?.x!;
      acl.y = event.accelerationIncludingGravity?.y!;
      console.log(JSON.stringify(event));
      const screenAngle = screen.orientation.angle ?? 0;
      // app.innerText = `aclX: ${acl.x}\naclY:${acl.y}\n`;
      if (screenAngle == 0) {
        scene.gravity = {
          x: -acl.x,
          y: -acl.y
        };
      } else if (screenAngle == 270) {
        scene.gravity = {
          x: -acl.y,
          y: acl.x,
        }
      } else if (screenAngle == 90) {
        scene.gravity = {
          x: acl.y,
          y: -acl.x,
        }
      }
    })
  }

  scene.setupScene({
    simWidth, simHeight, density: 1000
  });
  setupGUI(scene);

  function simulate() {
    scene.fluid.simulate(scene);
  }

  function update() {
    simulate();
    draw();
    requestAnimationFrame(update);
  }
  update();
}

main();