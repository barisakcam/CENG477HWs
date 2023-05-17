#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/

Vec3 translatePoint(Vec3 point, Translation tr) {
	point.x += tr.tx;
	point.y += tr.ty;
	point.z += tr.tz;

	return point;
}

Vec3 scalePoint(Vec3 point, Scaling sc) {
	point.x *= sc.sx;
	point.y *= sc.sy;
	point.z *= sc.sz;

	return point;
}

Vec3 rotatePoint(Vec3 point, Rotation rt) {
	Vec3 v, w, u;
	Vec4 result;
	u = normalizeVec3(Vec3(rt.ux, rt.uy, rt.uz, -1));
	if (abs(u.x) > abs(u.y)) {
		if (abs(u.y) > abs(u.z)) {
			v.x = -u.y;
			v.y = u.x;
			v.z = 0; //smallest: z
		} else {
			v.x = -u.z;
			v.y = 0;
			v.z = u.x; //smallest: y
		}
	} else if (abs(u.x) > abs(u.z)) {
		v.x = -u.y;
		v.y = u.x;
		v.z = 0; //smallest: z
	} else {
		v.x = 0;
		v.y = -u.z;
		v.z = u.y; //smallest: x
	}
	v = normalizeVec3(v);
	w = crossProductVec3(u, v);
	w = normalizeVec3(w);

	Matrix4 M, M_I, Rx;

	M.val[0][0] = u.x;
	M.val[0][1] = u.y;
	M.val[0][2] = u.z;
	M.val[0][3] = 0;

	M.val[1][0] = v.x;
	M.val[1][1] = v.y;
	M.val[1][2] = v.z;
	M.val[1][3] = 0;

	M.val[2][0] = w.x;
	M.val[2][1] = w.y;
	M.val[2][2] = w.z;
	M.val[2][3] = 0;

	M.val[3][0] = 0;
	M.val[3][1] = 0;
	M.val[3][2] = 0;
	M.val[3][3] = 1;

	M_I.val[0][0] = u.x;
	M_I.val[0][1] = v.x;
	M_I.val[0][2] = w.x;
	M_I.val[0][3] = 0;

	M_I.val[1][0] = u.y;
	M_I.val[1][1] = v.y;
	M_I.val[1][2] = w.y;
	M_I.val[1][3] = 0;

	M_I.val[2][0] = u.z;
	M_I.val[2][1] = v.z;
	M_I.val[2][2] = w.z;
	M_I.val[2][3] = 0;

	M_I.val[3][0] = 0;
	M_I.val[3][1] = 0;
	M_I.val[3][2] = 0;
	M_I.val[3][3] = 1;

	Rx.val[0][0] = 1;
	Rx.val[0][1] = 0;
	Rx.val[0][2] = 0;
	Rx.val[0][3] = 0;

	Rx.val[1][0] = 0;
	Rx.val[1][1] = cos((rt.angle * M_PI) / 180);
	Rx.val[1][2] = -sin((rt.angle * M_PI) / 180);
	Rx.val[1][3] = 0;

	Rx.val[2][0] = 0;
	Rx.val[2][1] = sin((rt.angle * M_PI) / 180);
	Rx.val[2][2] = cos((rt.angle * M_PI) / 180);
	Rx.val[2][3] = 0;

	Rx.val[3][0] = 0;
	Rx.val[3][1] = 0;
	Rx.val[3][2] = 0;
	Rx.val[3][3] = 1;

	result = multiplyMatrixWithVec4(multiplyMatrixWithMatrix(multiplyMatrixWithMatrix(M_I, Rx), M), Vec4(point.x, point.y, point.z, 1, -1));
	
	return Vec3(result.x, result.y, result.z, point.colorId);
}

Vec3 cameraTransform(Vec3 point,Camera c) {
	Matrix4 Mcam;
	Vec4 result;

	Mcam.val[0][0] = c.u.x;
	Mcam.val[0][1] = c.u.y;
	Mcam.val[0][2] = c.u.z;
	Mcam.val[0][3] = -(c.u.x * c.pos.x + c.u.y * c.pos.y + c.u.z * c.pos.z);

	Mcam.val[1][0] = c.v.x;
	Mcam.val[1][1] = c.v.y;
	Mcam.val[1][2] = c.v.z;
	Mcam.val[1][3] = -(c.v.x * c.pos.x + c.v.y * c.pos.y + c.v.z * c.pos.z);;

	Mcam.val[2][0] = c.w.x;
	Mcam.val[2][1] = c.w.y;
	Mcam.val[2][2] = c.w.z;
	Mcam.val[2][3] = -(c.w.x * c.pos.x + c.w.y * c.pos.y + c.w.z * c.pos.z);;

	Mcam.val[3][0] = 0;
	Mcam.val[3][1] = 0;
	Mcam.val[3][2] = 0;
	Mcam.val[3][3] = 1;

	result = multiplyMatrixWithVec4(Mcam, Vec4(point.x, point.y, point.z, 1, -1));

	return Vec3(result.x, result.y, result.z, point.colorId);
}

Vec4 projection(Vec3 point, Camera c) {
	Matrix4 Morth, Mper;
	Vec4 result;

	if (c.projectionType == 0) {
		Morth.val[0][0] = 2 / (c.right - c.left);
		Morth.val[0][1] = 0;
		Morth.val[0][2] = 0;
		Morth.val[0][3] = - (c.right + c.left) / (c.right - c.left);

		Morth.val[1][0] = 0;
		Morth.val[1][1] = 2 / (c.top - c.bottom);
		Morth.val[1][2] = 0;
		Morth.val[1][3] = - (c.top + c. bottom) / (c.top - c.bottom);

		Morth.val[2][0] = 0;
		Morth.val[2][1] = 0;
		Morth.val[2][2] = -2 / (c.far - c.near);
		Morth.val[2][3] = - (c.far + c.near) / (c.far - c.near);

		Morth.val[3][0] = 0;
		Morth.val[3][1] = 0;
		Morth.val[3][2] = 0;
		Morth.val[3][3] = 1;

		result = multiplyMatrixWithVec4(Morth, Vec4(point.x, point.y, point.z, 1, -1));
	} else {
		Mper.val[0][0] = (2 * c.near) / (c.right - c.left);
		Mper.val[0][1] = 0;
		Mper.val[0][2] = (c.right + c.left) / (c.right - c.left);
		Mper.val[0][3] = 0;

		Mper.val[1][0] = 0;
		Mper.val[1][1] = (2 * c.near) / (c.top - c.bottom);;
		Mper.val[1][2] = (c.top + c.bottom) / (c.top - c.bottom);;
		Mper.val[1][3] = 0;

		Mper.val[2][0] = 0;
		Mper.val[2][1] = 0;
		Mper.val[2][2] = - (c.far + c.near) / (c.far - c.near);
		Mper.val[2][3] = - (2 * c.far * c.near) / (c.far - c.near);

		Mper.val[3][0] = 0;
		Mper.val[3][1] = 0;
		Mper.val[3][2] = -1;
		Mper.val[3][3] = 0;

		result = multiplyMatrixWithVec4(Mper, Vec4(point.x, point.y, point.z, 1, -1));
	}

	result.colorId = point.colorId;
	return result;
}

Vec3 projectionDivide(Vec4 v) {
	return Vec3(v.x / v.t, v.y / v.t, v.z / v.t, v.colorId);
}

Vec3 viewport(Vec3 point, Camera c) {
	Matrix4 Mvp;
	Vec4 result;

	Mvp.val[0][0] = c.horRes / 2;
	Mvp.val[0][1] = 0;
	Mvp.val[0][2] = 0;
	Mvp.val[0][3] = (c.horRes - 1) / 2;

	Mvp.val[1][0] = 0;
	Mvp.val[1][1] = c.verRes / 2;
	Mvp.val[1][2] = 0;
	Mvp.val[1][3] = (c.verRes - 1) / 2;

	Mvp.val[2][0] = 0;
	Mvp.val[2][1] = 0;
	Mvp.val[2][2] = 0.5;
	Mvp.val[2][3] = 0.5;

	Mvp.val[3][0] = 0;
	Mvp.val[3][1] = 0;
	Mvp.val[3][2] = 0;
	Mvp.val[3][3] = 1;

	result = multiplyMatrixWithVec4(Mvp, Vec4(point.x, point.y, point.z, 1, -1));

	return Vec3(result.x / result.t, result.y / result.t, result.z / result.t, point.colorId);
}

bool visible(double den, double num, double &t_e, double &t_l) {
	double t;
	if (den > 0) {
		t = num / den;
		if (t > t_l) return false;
		if (t > t_e) t_e = t;
	} else if (den < 0) {
		t = num / den;
		if (t < t_e) return false;
		if (t < t_l) t_l = t;
	} else if (num > 0) return false;
	return true;
}

bool clippingLiangBarsky(Vec4 &point0, Vec4 &point1, Color &c0, Color &c1, Camera c) {
	double t_e = 0, t_l = 1;
	double d_x = point1.x - point0.x, d_y = point1.y - point0.y, d_z = point1.z - point0.z;
	double minw, maxw;
	Color dc;

	if (abs(point0.t) > abs(point1.t)) {
		maxw = abs(point0.t);
		minw = -abs(point0.t);
	} else {
		maxw = abs(point1.t);
		minw = -abs(point1.t);
	}

	dc.r = (c1.r - c0.r);
	dc.g = (c1.g - c0.g);
	dc.b = (c1.b - c0.b);

	bool result = false;
	if (visible(d_x, minw - point0.x, t_e, t_l)) {
		if (visible(-d_x, point0.x - maxw, t_e, t_l)) {
			if (visible(d_y, minw - point0.y, t_e, t_l)) {
				if (visible(-d_y, point0.y - maxw, t_e, t_l)) {
					if (visible(d_z, minw - point0.z, t_e, t_l)) {
						if (visible(-d_z, point0.z - maxw, t_e, t_l)) {
							result = true;
							if (t_l < 1) {
								point1.x = point0.x + d_x * t_l;
								point1.y = point0.y + d_y * t_l;
								point1.z = point0.z + d_z * t_l;
								c1.r = c0.r + dc.r * t_l;
								c1.g = c0.g + dc.g * t_l;
								c1.b = c0.b + dc.b * t_l;
							}
							if (t_e > 0) {
								point0.x = point0.x + d_x * t_e;
								point0.y = point0.y + d_y * t_e;
								point0.z = point0.z + d_z * t_e;
								c0.r = c0.r + dc.r * t_e;
								c0.g = c0.g + dc.g * t_e;
								c0.b = c0.b + dc.b * t_e;
							}
						}
					}
				}
			}
		}
	}
	return result;
}

bool backfaceCulling(Vec3 p1, Vec3 p2, Vec3 p3, Camera c) {
	Vec3 A, B, N, V;
	A = subtractVec3(p2, p1);
	B = subtractVec3(p3, p1);

	N.x = A.y * B.z - A.z * B.y;
	N.y = A.z * B.x - A.x * B.z;
	N.z = A.x * B.y - A.y * B.x;

	V.x = (p1.x + p2.x + p3.x) / 3;
	V.y = (p1.y + p2.y + p3.y) / 3;
	V.z = (p1.z + p2.z + p3.z) / 3;

	if (dotProductVec3(N, V) < 0) {
		return false;
	}
	return true;
}

void drawLow(int x0, int y0, int x1, int y1, Color c0, Color c1, vector< vector<Color> > &image) {
	int dx, dy, yi, D, y;
	Color c, dc;

	dx = x1 - x0;
	dy = y1 - y0;
	yi = 1;
	c = c0;
	dc.r = (c1.r - c0.r) / (x1 - x0);
	dc.g = (c1.g - c0.g) / (x1 - x0);
	dc.b = (c1.b - c0.b) / (x1 - x0);

	if (dy < 0) {
		yi = -1;
		dy = -dy;
	}
	D = (2 * dy) - dx;
	y = y0;

	for (int x = x0; x <= x1; x ++) {
		if (x >= 0 && x < image.size() && y >= 0 && y < image[0].size()) image[x][y] = c;
		if (D > 0) {
			y += yi;
			D += 2 * (dy - dx);
		} else {
			D += 2 * dy;
		}
		c.r += dc.r;
		c.g += dc.g;
		c.b += dc.b;
	}
}

void drawHigh(int x0, int y0, int x1, int y1, Color c0, Color c1, vector< vector<Color> > &image) {
	int dx, dy, xi, D, x;
	Color c, dc;

	dx = x1 - x0;
	dy = y1 - y0;
	xi = 1;
	c = c0;
	dc.r = (c1.r - c0.r) / (y1 - y0);
	dc.g = (c1.g - c0.g) / (y1 - y0);
	dc.b = (c1.b - c0.b) / (y1 - y0);

	if (dx < 0) {
		xi = -1;
		dx = -dx;
	}
	D = (2 * dx) - dy;
	x = x0;

	for (int y = y0; y <= y1; y ++) {
		if (x >= 0 && x < image.size() && y >= 0 && y < image[0].size()) image[x][y] = c;
		if (D > 0) {
			x += xi;
			D += 2 * (dx - dy);
		} else {
			D += 2 * dx;
		}
		c.r += dc.r;
		c.g += dc.g;
		c.b += dc.b;
	}
}

void lineRasterization(Vec3 p0, Vec3 p1, Color c0, Color c1, vector< vector<Color> > &image) {
	int p0_x, p0_y, p1_x, p1_y;

	Color c, dc;

	p0_x = (int) p0.x;
	p0_y = (int) p0.y;
	p1_x = (int) p1.x;
	p1_y = (int) p1.y;
	
	if (abs(p1_y - p0_y) < abs(p1_x - p0_x)) {
		if (p0_x > p1_x) {
			drawLow(p1_x, p1_y, p0_x, p0_y, c1, c0, image);
		} else {
			drawLow(p0_x, p0_y, p1_x, p1_y, c0, c1, image);
		}
	} else {
		if (p0_y > p1_y) {
			drawHigh(p1_x, p1_y, p0_x, p0_y, c1, c0, image);
		} else {
			drawHigh(p0_x, p0_y, p1_x, p1_y, c0, c1, image);
		}
	}
}

int f(int x, int y, int x_0, int y_0, int x_1, int y_1) {
	return x * (y_0 - y_1) + y * (x_1 - x_0) + x_0 * y_1 - y_0 * x_1;
}

void triangleRasterization(Vec3 p0, Vec3 p1, Vec3 p2, Color c0, Color c1, Color c2, vector< vector<Color> > &image) {
	int p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, y_min, y_max, x_min, x_max;
	double alph, beta, gama;
	Color c;
	p0_x = (int) p0.x;
	p0_y = (int) p0.y;
	p1_x = (int) p1.x;
	p1_y = (int) p1.y;
	p2_x = (int) p2.x;
	p2_y = (int) p2.y;

	y_min = p0_y;
	if (p1_y < y_min) y_min = p1_y;
	if (p2_y < y_min) y_min = p2_y;
	y_max = p0_y;
	if (p1_y > y_max) y_max = p1_y;
	if (p2_y > y_max) y_max = p2_y;
	x_min = p0_x;
	if (p1_x < x_min) x_min = p1_x;
	if (p2_x < x_min) x_min = p2_x;
	x_max = p0_x;
	if (p1_x > x_max) x_max = p1_x;
	if (p2_x > x_max) x_max = p2_x;

	for (int y = y_min; y <= y_max; y ++) {
		for (int x = x_min; x <= x_max; x ++) {
			alph = (double) f(x, y, p1_x, p1_y, p2_x, p2_y) / f(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y);
			beta = (double) f(x, y, p2_x, p2_y, p0_x, p0_y) / f(p1_x, p1_y, p2_x, p2_y, p0_x, p0_y);
			gama = (double) f(x, y, p0_x, p0_y, p1_x, p1_y) / f(p2_x, p2_y, p0_x, p0_y, p1_x, p1_y);
			if (alph >= 0 && beta >= 0 && gama >= 0) {
				c.r = alph * c0.r + beta * c1.r + gama * c2.r;
				c.g = alph * c0.g + beta * c1.g + gama * c2.g;
				c.b = alph * c0.b + beta * c1.b + gama * c2.b;
				if (x >= 0 && x < image.size() && y >= 0 && y < image[0].size()) image[x][y] = c;
			}
		}
	}
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	for (int i = 0; i < meshes.size(); i ++) {
		for (int j = 0; j < meshes[i]->triangles.size(); j ++) {
			Vec3 p0, p1, p2;
			Vec4 p0_w, p1_w, p2_w;
			p0 = *vertices[meshes[i]->triangles[j].getFirstVertexId() - 1];
			p1 = *vertices[meshes[i]->triangles[j].getSecondVertexId() - 1];
			p2 = *vertices[meshes[i]->triangles[j].getThirdVertexId() - 1];

			for (int k = 0; k < meshes[i]->numberOfTransformations; k ++) {
				if (meshes[i]->transformationTypes[k] == 't') {
					p0 = translatePoint(p0, *translations[meshes[i]->transformationIds[k] - 1]);
					p1 = translatePoint(p1, *translations[meshes[i]->transformationIds[k] - 1]);
					p2 = translatePoint(p2, *translations[meshes[i]->transformationIds[k] - 1]);
				} else if (meshes[i]->transformationTypes[k] == 's') {
					p0 = scalePoint(p0, *scalings[meshes[i]->transformationIds[k] - 1]);
					p1 = scalePoint(p1, *scalings[meshes[i]->transformationIds[k] - 1]);
					p2 = scalePoint(p2, *scalings[meshes[i]->transformationIds[k] - 1]);
				} else if (meshes[i]->transformationTypes[k] == 'r') {
					p0 = rotatePoint(p0, *rotations[meshes[i]->transformationIds[k] - 1]);
					p1 = rotatePoint(p1, *rotations[meshes[i]->transformationIds[k] - 1]);
					p2 = rotatePoint(p2, *rotations[meshes[i]->transformationIds[k] - 1]);
				}
			}
			
			p0 = cameraTransform(p0, *camera);
			p1 = cameraTransform(p1, *camera);
			p2 = cameraTransform(p2, *camera);
			p0_w = projection(p0, *camera);
			p1_w = projection(p1, *camera);
			p2_w = projection(p2, *camera);

			if (cullingEnabled && backfaceCulling(p0, p1, p2, *camera)) continue;

			if (meshes[i]->type == 0) {
				Vec4 l01_0 = p0_w, l01_1 = p1_w, l12_0 = p1_w, l12_1 = p2_w, l02_0 = p0_w, l02_1 = p2_w;
				Color l01_0_c = *colorsOfVertices[p0.colorId - 1], l01_1_c = *colorsOfVertices[p1.colorId - 1];
				Color l12_0_c = *colorsOfVertices[p1.colorId - 1], l12_1_c = *colorsOfVertices[p2.colorId - 1];
				Color l02_0_c = *colorsOfVertices[p0.colorId - 1], l02_1_c = *colorsOfVertices[p2.colorId - 1];
				if (clippingLiangBarsky(l01_0, l01_1, l01_0_c, l01_1_c, *camera)) {
					Vec3 l0, l1;
					l0 = projectionDivide(l01_0);
					l1 = projectionDivide(l01_1);
					l0 = viewport(l0, *camera);
					l1 = viewport(l1, *camera);
					lineRasterization(l0, l1, l01_0_c, l01_1_c, image);
				}
				if (clippingLiangBarsky(l12_0, l12_1, l12_0_c, l12_1_c, *camera)) {
					Vec3 l0, l1;
					l0 = projectionDivide(l12_0);
					l1 = projectionDivide(l12_1);
					l0 = viewport(l0, *camera);
					l1 = viewport(l1, *camera);
					lineRasterization(l0, l1, l12_0_c, l12_1_c, image);
				}
				if (clippingLiangBarsky(l02_0, l02_1, l02_0_c, l02_1_c, *camera)) {
					Vec3 l0, l1;
					l0 = projectionDivide(l02_0);
					l1 = projectionDivide(l02_1);
					l0 = viewport(l0, *camera);
					l1 = viewport(l1, *camera);
					lineRasterization(l0, l1, l02_0_c, l02_1_c, image);
				}
			} else {
				Vec3 l0, l1, l2;
				l0 = projectionDivide(p0_w);
				l1 = projectionDivide(p1_w);
				l2 = projectionDivide(p2_w);
				l0 = viewport(l0, *camera);
				l1 = viewport(l1, *camera);
				l2 = viewport(l2, *camera);
				triangleRasterization(l0, l1, l2, *colorsOfVertices[p0.colorId - 1], *colorsOfVertices[p1.colorId - 1], *colorsOfVertices[p2.colorId - 1], image);
			}
		}
	}
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL) {
		str = pElement->GetText();
		
		if (strcmp(str, "enabled") == 0) {
			cullingEnabled = true;
		}
		else {
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0) {
			cam->projectionType = 0;
		}
		else {
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0) {
			mesh->type = 0;
		}
		else {
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
        str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);
			
			if (result != EOF) {
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}