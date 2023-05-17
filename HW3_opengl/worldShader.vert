#version 420

layout (location = 0) in vec3 VertexPosition;
layout (location = 1) in vec3 VertexNormal;
layout (location = 2) in vec2 VertexTex;

uniform vec3 lightPosition;
uniform vec3 cameraPosition;

uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;
uniform mat4 NormalMatrix;
uniform mat4 MVP;

uniform sampler2D TexColor;
uniform sampler2D TexGrey;
uniform float textureOffset;

uniform float heightFactor;
uniform float imageWidth;
uniform float imageHeight;

out Data
{
    vec3 Position;
    vec3 Normal;
    vec2 TexCoord;
} data;


out vec3 LightVector;// Vector from Vertex to Light;
out vec3 CameraVector;// Vector from Vertex to Camera;

void main()
{
    vec2 textureCoordinate = vec2(VertexTex.x, VertexTex.y);
    vec4 texColor = texture(TexGrey, textureCoordinate);
    float height = heightFactor * texColor.x;
    vec4 pos, normal;

    normal = MVP * vec4(VertexNormal, 1.0);
    pos = MVP * vec4(VertexPosition.x, VertexPosition.y, VertexPosition.z, 1.0);
    pos.x += normal.x * height;
    pos.y += normal.y * height;
    pos.z += normal.z * height;
    pos = ProjectionMatrix * ViewMatrix * pos;

    data.TexCoord = VertexTex;
    data.Position.x = pos.x;
    data.Position.y = pos.y;
    data.Position.z = pos.z;
    data.Normal.x = normal.x;
    data.Normal.y = normal.y;
    data.Normal.z = normal.z;

    LightVector = normalize(lightPosition - data.Position);
    CameraVector = normalize(cameraPosition - data.Position);
    
    gl_Position = pos;
}