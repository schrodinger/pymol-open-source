uniform sampler2D diffuse;
varying vec2 texcoords;

void main() {
  gl_FragColor = texture2D(diffuse, texcoords);
}
