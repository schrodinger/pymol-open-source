attribute vec3 position;
attribute vec2 texcoords_in;
varying vec2 texcoords;

void main() {
  texcoords = texcoords_in;
  gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 1);
}
