#include <math.h>
#include "parser.h"
#include "ppm.h"

typedef parser::Vec3f Vec;

int nx, ny;
double dist, left, right, bottom, top;
Vec e, gaze, v, w, u;
parser::Scene scene;

typedef struct {
    Vec o, d;
} Ray;

double max(double a, double b) {
    if (a > b) return a;
    else return b;
}

double min(double a, double b) {
    if (a < b) return a;
    else return b;
}

double length(Vec a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

Vec scalar_mult(Vec a, double s) {
    Vec result;
    result.x = a.x * s;
    result.y = a.y * s;
    result.z = a.z * s;
    return result;
}

Vec add(Vec a, Vec b) {
    Vec result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
	return result;
}

double dot_product(Vec a, Vec b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec cross_product(Vec a, Vec b) {
    Vec result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

Vec normalize(Vec a) {
    return scalar_mult(a, 1.0 / sqrt(dot_product(a, a)));
}

Ray create_ray(int i, int j) {
    Ray result;
    double su, sv;
    Vec m, q, s;
    
    su = (i + 0.5) * (right - left) / nx;
    sv = (j + 0.5) * (top -bottom) / ny;

    m = add(e, scalar_mult(gaze, dist));
    q = add(m, add(scalar_mult(u, left), scalar_mult(v, top)));
    s = add(q, add(scalar_mult(u, su), scalar_mult(v, -sv)));

    result.o = e;
    result.d = add(s, scalar_mult(e, -1));

    return result;
}

double determinant(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g); 
}

Vec tri_normal(parser::Face a) {
    Vec result;
    Vec v0, v1, v2;
    v0 = scene.vertex_data[a.v0_id - 1];
    v1 = scene.vertex_data[a.v1_id - 1];
    v2 = scene.vertex_data[a.v2_id - 1];
    result.x = ((v1.y - v0.y) * (v2.z - v0.z)) - ((v1.z - v0.z) * (v2.y - v0.y));
    result.y = ((v1.z - v0.z) * (v2.x - v0.x)) - ((v1.x - v0.x) * (v2.z - v0.z));
    result.z = ((v1.x - v0.x) * (v2.y - v0.y)) - ((v1.y - v0.y) * (v2.x - v0.x));
    return result;
}

double intersect_triangle(Ray r, Vec a, Vec b, Vec c) {
    double detA = determinant(a.x - b.x, a.y - b.y, a.z - b.z, a.x - c.x, a.y - c.y, a.z - c.z, r.d.x, r.d.y, r.d.z);
    double detB = determinant(a.x - r.o.x, a.y - r.o.y, a.z - r.o.z, a.x - c.x, a.y - c.y, a.z - c.z, r.d.x, r.d.y, r.d.z);
    double detD = determinant(a.x - b.x, a.y - b.y, a.z - b.z, a.x - r.o.x, a.y - r.o.y, a.z - r.o.z, r.d.x, r.d.y, r.d.z);
    double detT = determinant(a.x - b.x, a.y - b.y, a.z - b.z, a.x - c.x, a.y - c.y, a.z - c.z, a.x - r.o.x, a.y - r.o.y, a.z - r.o.z);
    double beta = detB / detA;
    double delta = detD / detA;
    double t = detT / detA;

    double minT = -1, maxT = 99999999;

    if (minT <= t && maxT >= t && beta + delta <= 1 && 0 <= beta && 0 <= delta) {
        return t;
    } else {
        return -1;
    }
}

double intersect_sphere(Ray r, parser::Sphere s) {
    double A, B, C;
    double delta;
    Vec c;
    c = scene.vertex_data[s.center_vertex_id - 1];

    double t, t1, t2;

    C = (r.o.x - c.x) * (r.o.x - c.x) + (r.o.y - c.y) * (r.o.y - c.y) + (r.o.z - c.z) * (r.o.z - c.z) - s.radius * s.radius;
    B = 2 * r.d.x * (r.o.x - c.x) + 2 * r.d.y * (r.o.y - c.y) + 2 * r.d.z * (r.o.z - c.z);
    A = r.d.x * r.d.x + r.d.y * r.d.y + r.d.z * r.d.z;
    delta = B * B - 4 * A * C;

    if (delta < 0) return -1;
    else if (delta == 0) {
        t = -B / (2 * A);
    } else {
        t1 = (-B + sqrt(delta)) / (2 * A);
        t2 = (-B - sqrt(delta)) / (2 * A);
        if (t1 < t2) t = t1;
        else t = t2;
    }

    return t;
}

bool compute_shadow(Ray r, double dist) {
    double t;

    for (int i = 0; i < scene.spheres.size(); i ++) {
        t = intersect_sphere(r, scene.spheres[i]);
        if (t > 0 && t < dist) return true;
    }

    for (int i = 0; i < scene.triangles.size(); i++) {
        t = intersect_triangle(r, scene.vertex_data[scene.triangles[i].indices.v0_id - 1], scene.vertex_data[scene.triangles[i].indices.v1_id - 1], scene.vertex_data[scene.triangles[i].indices.v2_id - 1]);
        if (t > 0 && t < dist) return true;
    }

    for (int i = 0; i < scene.meshes.size(); i ++) {
        for (int j = 0; j < scene.meshes[i].faces.size(); j ++) {
            t = intersect_triangle(r, scene.vertex_data[scene.meshes[i].faces[j].v0_id - 1], scene.vertex_data[scene.meshes[i].faces[j].v1_id - 1], scene.vertex_data[scene.meshes[i].faces[j].v2_id - 1]);
            if (t > 0 && t < dist) return true;
        }
    }

    return false;
}

Vec compute_color(Ray r, int recursion) {
    Ray shadow_ray;
    Vec d, s, m1, m2;
    Vec I;
    Vec I_pos;
    Vec c;
    Vec n;
    Vec w_i, w_i_norm, h, w_r;
    Vec L_d, L_s, L_m;
    double I_dist;
    double minT = 99999999;
    double t;
    int minI = -1, minJ = -1;
    int hit = 0;

    c.x = 0;
    c.y = 0;
    c.z = 0;

    if (recursion < 0) return c;

    for (int i = 0; i < scene.spheres.size(); i ++) {
        t = intersect_sphere(r, scene.spheres[i]);
        if (t < minT && t >= 0) {
            minI = i;
            minT = t;
            hit = 1;
        }
    }

    for (int i = 0; i < scene.triangles.size(); i++) {
        t = intersect_triangle(r, scene.vertex_data[scene.triangles[i].indices.v0_id - 1], scene.vertex_data[scene.triangles[i].indices.v1_id - 1], scene.vertex_data[scene.triangles[i].indices.v2_id - 1]);
        if (t < minT && t >= 0) {
            minI = i;
            minT = t;
            hit = 2;
        }
    }

    for (int i = 0; i < scene.meshes.size(); i ++) {
        for (int j = 0; j < scene.meshes[i].faces.size(); j ++) {
            t = intersect_triangle(r, scene.vertex_data[scene.meshes[i].faces[j].v0_id - 1], scene.vertex_data[scene.meshes[i].faces[j].v1_id - 1], scene.vertex_data[scene.meshes[i].faces[j].v2_id - 1]);
            if (t < minT && t >= 0) {
                minI = i;
                minJ = j;
                minT = t;
                hit = 3;
            }
        }
    }

    switch (hit) {
        case 0:
            if (recursion == scene.max_recursion_depth) {
                c.x = scene.background_color.x;
                c.y = scene.background_color.y;
                c.z = scene.background_color.z;
            }
            break;
        case 1:
            c.x = scene.materials[scene.spheres[minI].material_id - 1].ambient.x * scene.ambient_light.x;
            c.y = scene.materials[scene.spheres[minI].material_id - 1].ambient.y * scene.ambient_light.y;
            c.z = scene.materials[scene.spheres[minI].material_id - 1].ambient.z * scene.ambient_light.z;
            d = scene.materials[scene.spheres[minI].material_id - 1].diffuse;
            s = scene.materials[scene.spheres[minI].material_id - 1].specular;
            n = add(add(r.o, scalar_mult(r.d, minT)), scalar_mult(scene.vertex_data[scene.spheres[minI].center_vertex_id - 1], -1));
            n = normalize(n);

            for (int i = 0; i < scene.point_lights.size(); i ++) {
                I_pos = scene.point_lights[i].position;
                I = scene.point_lights[i].intensity;
                w_i = add(scalar_mult(add(r.o, scalar_mult(r.d, minT)), -1), I_pos);
                I_dist = length(w_i);
                w_i_norm = normalize(w_i);
                shadow_ray.o = add(add(r.o, scalar_mult(r.d, minT)), scalar_mult(n, scene.shadow_ray_epsilon));
                shadow_ray.d = w_i_norm;
                if (compute_shadow(shadow_ray, length(w_i))) continue;

                L_d.x = d.x * (I.x / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                L_d.y = d.y * (I.y / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                L_d.z = d.z * (I.z / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                
                h = add(w_i_norm, normalize(scalar_mult(r.d, -1)));
                h = normalize(h);
                if (dot_product(n, w_i_norm) > 0) {
                    L_s.x = s.x * (I.x / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.spheres[minI].material_id - 1].phong_exponent);
                    L_s.y = s.y * (I.y / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.spheres[minI].material_id - 1].phong_exponent);
                    L_s.z = s.z * (I.z / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.spheres[minI].material_id - 1].phong_exponent);
                } else {
                    L_s.x = 0;
                    L_s.y = 0;
                    L_s.z = 0;
                }

                c = add(c, L_d);
                c = add(c, L_s);
            }

            if (scene.materials[scene.spheres[minI].material_id - 1].is_mirror) {
                Ray mirror_ray;
                w_r = add(normalize(r.d), scalar_mult(n, 2 * dot_product(n, normalize(scalar_mult(r.d, -1)))));
                mirror_ray.o = add(add(r.o, scalar_mult(r.d, minT)), scalar_mult(n, scene.shadow_ray_epsilon));
                mirror_ray.d = w_r;
                m1 = scene.materials[scene.spheres[minI].material_id - 1].mirror;
                m2 = compute_color(mirror_ray, recursion - 1);
                L_m.x = m1.x * m2.x;
                L_m.y = m1.y * m2.y;
                L_m.z = m1.z * m2.z; 
                c = add(c, L_m);
            }
            break;
        case 2:
            c.x = scene.materials[scene.triangles[minI].material_id - 1].ambient.x * scene.ambient_light.x;
            c.y = scene.materials[scene.triangles[minI].material_id - 1].ambient.y * scene.ambient_light.y;
            c.z = scene.materials[scene.triangles[minI].material_id - 1].ambient.z * scene.ambient_light.z;
            d = scene.materials[scene.triangles[minI].material_id - 1].diffuse;
            s = scene.materials[scene.triangles[minI].material_id - 1].specular;
            n = tri_normal(scene.triangles[minI].indices);
            n = normalize(n);

            for (int i = 0; i < scene.point_lights.size(); i ++) {
                I_pos = scene.point_lights[i].position;
                I = scene.point_lights[i].intensity;
                w_i = add(scalar_mult(add(r.o, scalar_mult(r.d, minT)), -1), I_pos);
                I_dist = length(w_i);
                w_i_norm = normalize(w_i);
                shadow_ray.o = add(add(r.o, scalar_mult(r.d, minT)), scalar_mult(n, scene.shadow_ray_epsilon));
                shadow_ray.d = w_i_norm;
                if (compute_shadow(shadow_ray, length(w_i))) continue;

                L_d.x = d.x * (I.x / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                L_d.y = d.y * (I.y / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                L_d.z = d.z * (I.z / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                
                h = add(w_i_norm, normalize(scalar_mult(r.d, -1)));
                h = normalize(h);
                if (dot_product(n, w_i_norm) > 0) {
                    L_s.x = s.x * (I.x / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.triangles[minI].material_id - 1].phong_exponent);
                    L_s.y = s.y * (I.y / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.triangles[minI].material_id - 1].phong_exponent);
                    L_s.z = s.z * (I.z / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.triangles[minI].material_id - 1].phong_exponent);
                } else {
                    L_s.x = 0;
                    L_s.y = 0;
                    L_s.z = 0;
                }

                c = add(c, L_d);
                c = add(c, L_s);
            }

            if (scene.materials[scene.triangles[minI].material_id - 1].is_mirror) {
                Ray mirror_ray;
                w_r = add(normalize(r.d), scalar_mult(n, 2 * dot_product(n, normalize(scalar_mult(r.d, -1)))));
                mirror_ray.o = add(add(r.o, scalar_mult(r.d, minT)), scalar_mult(n, scene.shadow_ray_epsilon));
                mirror_ray.d = w_r;
                m1 = scene.materials[scene.triangles[minI].material_id - 1].mirror;
                m2 = compute_color(mirror_ray, recursion - 1);
                L_m.x = m1.x * m2.x;
                L_m.y = m1.y * m2.y;
                L_m.z = m1.z * m2.z; 
                c = add(c, L_m);
            }
            break;
        case 3:
            c.x = scene.materials[scene.meshes[minI].material_id - 1].ambient.x * scene.ambient_light.x;
            c.y = scene.materials[scene.meshes[minI].material_id - 1].ambient.y * scene.ambient_light.y;
            c.z = scene.materials[scene.meshes[minI].material_id - 1].ambient.z * scene.ambient_light.z;
            d = scene.materials[scene.meshes[minI].material_id - 1].diffuse;
            s = scene.materials[scene.meshes[minI].material_id - 1].specular;
            n = tri_normal(scene.meshes[minI].faces[minJ]);
            n = normalize(n);

            for (int i = 0; i < scene.point_lights.size(); i ++) {
                I_pos = scene.point_lights[i].position;
                I = scene.point_lights[i].intensity;
                w_i = add(scalar_mult(add(r.o, scalar_mult(r.d, minT)), -1), I_pos);
                I_dist = length(w_i);
                w_i_norm = normalize(w_i);
                shadow_ray.o = add(add(r.o, scalar_mult(r.d, minT)), scalar_mult(n, scene.shadow_ray_epsilon));
                shadow_ray.d = w_i_norm;
                if (compute_shadow(shadow_ray, length(w_i))) continue;

                L_d.x = d.x * (I.x / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                L_d.y = d.y * (I.y / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));
                L_d.z = d.z * (I.z / (I_dist * I_dist)) * max(0, dot_product(n, w_i_norm));

                h = add(w_i_norm, normalize(scalar_mult(r.d, -1)));
                h = normalize(h);
                if (dot_product(n, w_i_norm) > 0) {
                    L_s.x = s.x * (I.x / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.meshes[minI].material_id - 1].phong_exponent);
                    L_s.y = s.y * (I.y / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.meshes[minI].material_id - 1].phong_exponent);
                    L_s.z = s.z * (I.z / (I_dist * I_dist)) * pow(max(0, dot_product(n, h)), scene.materials[scene.meshes[minI].material_id - 1].phong_exponent);
                } else {
                    L_s.x = 0;
                    L_s.y = 0;
                    L_s.z = 0;
                }

                c = add(c, L_d);
                c = add(c, L_s);
            }

            if (scene.materials[scene.meshes[minI].material_id - 1].is_mirror) {
                Ray mirror_ray;
                w_r = add(normalize(r.d), scalar_mult(n, 2 * dot_product(n, normalize(scalar_mult(r.d, -1)))));
                mirror_ray.o = add(add(r.o, scalar_mult(r.d, minT)), scalar_mult(n, scene.shadow_ray_epsilon));
                mirror_ray.d = w_r;
                m1 = scene.materials[scene.meshes[minI].material_id - 1].mirror;
                m2 = compute_color(mirror_ray, recursion - 1);
                L_m.x = m1.x * m2.x;
                L_m.y = m1.y * m2.y;
                L_m.z = m1.z * m2.z;
                c = add(c, L_m);
            }
            break;
    }

    c.x = min(c.x, 255);
    c.y = min(c.y, 255);
    c.z = min(c.z, 255);

    return c;
}

int main(int argc, char* argv[])
{
    scene.loadFromXml(argv[1]);

    for (int camera_counter = 0; camera_counter < scene.cameras.size(); camera_counter ++) {
        nx = scene.cameras[camera_counter].image_width;
        ny = scene.cameras[camera_counter].image_height;
        dist = scene.cameras[camera_counter].near_distance;
        left = scene.cameras[camera_counter].near_plane.x;
        right = scene.cameras[camera_counter].near_plane.y;
        bottom = scene.cameras[camera_counter].near_plane.z;
        top = scene.cameras[camera_counter].near_plane.w;
        e = scene.cameras[camera_counter].position;
        gaze = scene.cameras[camera_counter].gaze;
        v = scene.cameras[camera_counter].up;
        w = scalar_mult(gaze, -1);
        u = cross_product(v, w);

        unsigned char * image = new unsigned char [nx * ny * 3];

        int k = 0;
        for (int j = 0; j < ny; j ++) {
            for (int i = 0; i < nx; i ++) {
                Ray new_ray = create_ray(i, j);

                Vec pixel = add(new_ray.o, new_ray.d);
                Vec color = compute_color(new_ray, scene.max_recursion_depth);

                image[k ++] = (int) (color.x);
                image[k ++] = (int) (color.y);
                image[k ++] = (int) (color.z);
            }
        }
        write_ppm(scene.cameras[camera_counter].image_name.c_str(), image, nx, ny);
    }
}
