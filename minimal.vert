// Vertex Shader – file "minimal.vert"

#version 130

in vec3 vertexPosition_modelspace;
uniform mat4 MVP;
in  vec3 in_Color;
out vec3 ex_Color;

void main(){
    vec4 v = vec4(vertexPosition_modelspace,1);
    gl_Position = MVP * v;
    ex_Color = in_Color;
}
