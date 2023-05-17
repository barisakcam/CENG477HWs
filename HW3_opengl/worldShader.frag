#version 420

in Data
{
    vec3 Position;
    vec3 Normal;
    vec2 TexCoord;
} data;
in vec3 LightVector;
in vec3 CameraVector;

uniform vec3 lightPosition;
uniform sampler2D TexColor;
uniform sampler2D MoonTexColor;
uniform sampler2D TexGrey;
uniform float textureOffset;

out vec4 FragColor;

vec3 ambientReflectenceCoefficient = vec3(0.5, 0.5, 0.5);
vec3 ambientLightColor = vec3(0.6, 0.6, 0.6);
vec3 specularReflectenceCoefficient= vec3(1.0f);
vec3 specularLightColor = vec3(1.0f);
float SpecularExponent = 10;
vec3 diffuseReflectenceCoefficient= vec3(1.0f);
vec3 diffuseLightColor = vec3(1.0f);


void main()
{
    vec2 textureCoordinate = vec2(data.TexCoord.x, data.TexCoord.y);
    vec4 texColor = texture(TexColor, textureCoordinate);
    vec3 halfVector = normalize(LightVector + CameraVector);

    vec3 ambient;
    ambient.x = ambientReflectenceCoefficient.x * ambientLightColor.x * texColor.x;
    ambient.y = ambientReflectenceCoefficient.y * ambientLightColor.y * texColor.y;
    ambient.z = ambientReflectenceCoefficient.z * ambientLightColor.z * texColor.z;

    vec3 diffuse;
    diffuse.x = texColor.x * max(0, dot(LightVector, data.Normal)) * diffuseLightColor.x;
    diffuse.y = texColor.y * max(0, dot(LightVector, data.Normal)) * diffuseLightColor.y;
    diffuse.z = texColor.z * max(0, dot(LightVector, data.Normal)) * diffuseLightColor.z;

    vec3 spec;
    spec.x = specularReflectenceCoefficient.x * pow(max(0, dot(data.Normal, halfVector)), SpecularExponent) * texColor.x;
    spec.y = specularReflectenceCoefficient.y * pow(max(0, dot(data.Normal, halfVector)), SpecularExponent) * texColor.y;
    spec.z = specularReflectenceCoefficient.z * pow(max(0, dot(data.Normal, halfVector)), SpecularExponent) * texColor.z;

    FragColor = vec4(ambient+diffuse+spec, 1.0f);
}
