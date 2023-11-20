
async function main() {
  const app = document.querySelector<HTMLDivElement>('#app')!;
  const canvas = document.querySelector<HTMLCanvasElement>('#container')!;
  const gl = canvas.getContext('webgl')!;
  const vertexShaderSource = `
        attribute vec2 a_position;
        uniform vec2 u_resolution;
        uniform vec2 u_translation;
        void main() {
          vec2 position = a_position + u_translation;
          vec2 clipSpace = (position / u_resolution) * 2.0 - 1.0;
          gl_Position = vec4(clipSpace * vec2(1, -1), 0, 1);
        }
      `;

  const fragmentShaderSource = `
        precision mediump float;
        uniform vec4 u_color;
        void main() {
          gl_FragColor = u_color;
        }
      `;

  function createShader(gl: WebGLRenderingContext, type: number, source: WebGLShader & string) {
    const shader = gl.createShader(type)!;
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      console.error('An error occurred compiling the shaders:', gl.getShaderInfoLog(shader));
      gl.deleteShader(shader);
      return null;
    }
    return shader;
  }

  const vertexShader = createShader(gl, gl.VERTEX_SHADER, vertexShaderSource)!;
  const fragmentShader = createShader(gl, gl.FRAGMENT_SHADER, fragmentShaderSource)!;

  const program = gl.createProgram()!;
  gl.attachShader(program, vertexShader);
  gl.attachShader(program, fragmentShader);
  gl.linkProgram(program);

  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    console.error('Unable to initialize the shader program:', gl.getProgramInfoLog(program));
  }

  gl.useProgram(program);

  const positionBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);

  const positionAttributeLocation = gl.getAttribLocation(program, 'a_position');
  gl.enableVertexAttribArray(positionAttributeLocation);
  gl.vertexAttribPointer(positionAttributeLocation, 2, gl.FLOAT, false, 0, 0);

  const resolutionUniformLocation = gl.getUniformLocation(program, 'u_resolution');
  const translationUniformLocation = gl.getUniformLocation(program, 'u_translation');
  const colorUniformLocation = gl.getUniformLocation(program, 'u_color');

  gl.uniform2f(resolutionUniformLocation, canvas.width, canvas.height);

  const circle = {
    x: canvas.width / 2,
    y: canvas.height / 2,
    radius: 30,
    vx: 2,
    vy: 6,
    color: [1, 0, 0, 1],
  };

  // draw a circle: using many segments
  const numSegments = 64;
  // calculate the position of every triangle vertexs
  const positions = [0, 0];
  for (let i = 0; i <= numSegments; i++) {
    const angle = (i / numSegments) * 2 * Math.PI;
    positions.push(Math.cos(angle) * circle.radius, Math.sin(angle) * circle.radius);
  }

  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(positions), gl.STATIC_DRAW);


  const gravity = 9.8 / 50;
  // acceleration
  // -> x up
  // |
  // y up
  //
  let acl = {
    x: 0,
    y: gravity,
  };
  window.addEventListener('devicemotion', (event: DeviceMotionEvent) => {
    acl.x = -event.accelerationIncludingGravity?.x! / 100;
    acl.y = event.accelerationIncludingGravity?.y! / 100;
    app.innerText = `aclX: ${acl.x}\naclY:${acl.y}`;
  })

  const friction = 0.9;
  gl.clearColor(0, 0, 0, 1);

  function animate() {
    gl.clear(gl.COLOR_BUFFER_BIT);

    circle.vx += acl.x;
    circle.vy += acl.y;
    // console.log('acc:', acl.accX, acl.accY);

    circle.x += circle.vx;
    circle.y += circle.vy;

    // bounce
    if (circle.x + circle.radius > canvas.width) {
      circle.x = canvas.width - circle.radius;
      circle.vx = -circle.vx;
    }
    if (circle.x - circle.radius < 0) {
      circle.x = circle.radius;
      circle.vx = -circle.vx;
    }

    // bounce
    if (circle.y + circle.radius > canvas.height) {
      circle.y = canvas.height - circle.radius;
      circle.vy = -circle.vy * friction;
    }
    if (circle.y - circle.radius < 0) {
      circle.y = circle.radius;
      circle.vy = -circle.vy * friction;
    }

    gl.uniform2f(translationUniformLocation, circle.x, circle.y);
    gl.uniform4fv(colorUniformLocation, circle.color);

    gl.drawArrays(gl.TRIANGLE_FAN, 0, numSegments + 2);

    requestAnimationFrame(animate);
  }

  animate();
}


main();
