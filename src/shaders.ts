export const pointVertexShader = `
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

export const pointFragmentShader = `
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

export const meshVertexShader = `
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

export const meshFragmentShader = `
  precision mediump float;
  varying vec3 fragColor;

  void main() {
    gl_FragColor = vec4(fragColor, 1.0);
  }
`;

export function createShader(gl: WebGLRenderingContext, vsSource: WebGLShader & string, fsSource: WebGLShader & string) {
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