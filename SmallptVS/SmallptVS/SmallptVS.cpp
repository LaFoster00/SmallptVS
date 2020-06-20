#include <math.h>
#include <random>
#include <stdio.h>
#include <omp.h>

#define M_PI 3.1415926535897932384626433832795
//#define M_PI 0.5
#define Samples 1000
#define Depth 5
#define MaxDepth 20

inline double erand() {
	return (double)rand()
		/ (double)RAND_MAX;
}

struct Vec {

	double x, y, z;

	Vec(double x_ = 0, double y_ = 0, double z_ = 0)
	{
		x = x_;
		y = y_;
		z = z_;
	}

	Vec operator+(const Vec& b) const
	{
		return Vec(x + b.x, y + b.y, z + b.z);
	}

	Vec operator-(const Vec& b) const
	{
		return Vec(x - b.x, y - b.y, z - b.z);
	}

	Vec operator*(double b) const
	{
		return Vec(x * b, y * b, z * b);
	}

	Vec mult(const Vec& b) const
	{
		return Vec(x * b.x, y * b.y, z * b.z);
	}

	Vec& norm()
	{
		return *this = *this * (1 / sqrt(x * x + y * y + z * z));
	}

	double dot(const Vec& b) const
	{
		return x * b.x + y * b.y + z * b.z;
	}

	// cross:
	Vec operator%(Vec& b)
	{
		return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
	}

};

struct Ray
{
	Vec origin, dir;
	Ray(Vec origin_, Vec dir_)
	{
		origin = origin_;
		dir = dir_;
	}
};

// material types, used in radiance() 
enum class Refl_t
{
	DIFF,
	SPEC,
	REFR
};

struct Sphere {
	double rad;       // radius 
	Vec pos, emi, col;      // position, emission, color 
	Refl_t reflT;      // reflection type (DIFFuse, SPECular, REFRactive) 

	Sphere(double rad_, Vec pos_, Vec emi_, Vec col_, Refl_t reflT_)
	{
		rad = rad_;
		pos = pos_;
		emi = emi_;
		col = col_;
		reflT = reflT_;
	}

	// returns distance, 0 if nohit 
	double intersect(const Ray& r) const
	{
		Vec offsetPosition = pos - r.origin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
		double t;
		double eps = 1e-4, b = offsetPosition.dot(r.dir);
		double det = b * b - offsetPosition.dot(offsetPosition) + rad * rad;

		if (det < 0)
		{
			return 0;
		}
		else
		{
			det = sqrt(det);
		}

		return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
	}
};

Sphere spheres[] = {//Scene: radius, position, emission, color, material 
  Sphere(1e5, Vec(1e5 + 1,40.8,81.6),	Vec(),			Vec(.75,.65,.25),	Refl_t::DIFF),//Left 
  Sphere(1e5, Vec(-1e5 + 99,40.8,81.6),	Vec(),			Vec(.25,.25,.75),	Refl_t::DIFF),//Right 
  Sphere(1e5, Vec(50,40.8, 1e5),		Vec(),			Vec(.75,.75,.75),	Refl_t::DIFF),//Back 
  Sphere(1e5, Vec(50,40.8,-1e5 + 170),	Vec(),			Vec(),				Refl_t::DIFF),//Front 
  Sphere(1e5, Vec(50, 1e5, 81.6),		Vec(),			Vec(.75,.75,.75),	Refl_t::DIFF),//Botom 
  Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),	Vec(),			Vec(.75,.75,.75),	Refl_t::DIFF),//Top 
  Sphere(16.5,Vec(27,16.5,47),			Vec(),			Vec(1,1,1) * .999,	Refl_t::SPEC),//Miror 
  Sphere(16.5,Vec(73,16.5,78),			Vec(),			Vec(1,1,1) * .999,	Refl_t::REFR),//Glass 
  Sphere(600, Vec(50,681.6 - .27,81.6),	Vec(12,12,12),	Vec(),				Refl_t::DIFF) //Light 
};

inline double clamp(double x)
{
	return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x)
{
	return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

void output_ppm(const char* fileName, int width, int height, Vec* img)
{
	FILE* file = fopen(fileName, "w");         // Write image to PPM file. 
	fprintf(file, "P3\n%d %d\n%d\n", width, height, 255);
	for (int i = 0; i < width * height; i++)
		fprintf(file, "%d %d %d ", toInt(img[i].x), toInt(img[i].y), toInt(img[i].z));
}

inline bool intersectWithScene(const Ray& r, double& t, int& id)
{
	double n = sizeof(spheres) / sizeof(Sphere);
	double d;
	double inf = t = 1e20;
	for (int i = n; i--;)
	{
		if ((d = spheres[i].intersect(r)) && d < t)
		{
			t = d;
			id = i;
		}
	}
	return t < inf;
}

Vec radiance(const Ray& ray, int depth) {
	double hitDistance; // distance to intersection 
	int id = 0; // id of intersected object 
	if (!intersectWithScene(ray, hitDistance, id))
	{
		return Vec(); // if miss, return black
	}

	const Sphere& obj = spheres[id]; // the hit object 

	Vec hitPos = ray.origin + ray.dir * hitDistance;
	Vec hitNormal = (hitPos - obj.pos).norm();
	Vec hitNormalNegative = hitNormal.dot(ray.dir) < 0 ? hitNormal : hitNormal * -1;
	Vec hitColor = obj.col;

	double maxReflectionStrength = hitColor.x > hitColor.y&& hitColor.x > hitColor.z ? hitColor.x : hitColor.y > hitColor.z ? hitColor.y : hitColor.z; // max reflection strength (brightest color component)

	if (++depth > Depth)
	{
		if (erand() < maxReflectionStrength && depth <= MaxDepth)
		{
			hitColor = hitColor * (1 / maxReflectionStrength);
		}
		else
		{
			return obj.emi;
		}
	}//R.R. 

	if (obj.reflT == Refl_t::DIFF) // Ideal DIFFUSE reflection 
	{
		double r1 = 2 * M_PI * erand();
		double r2 = erand();
		double r2s = sqrt(r2);
		Vec w = hitNormalNegative;
		Vec u = ((abs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w % u;
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

		return obj.emi + hitColor.mult(radiance(Ray(hitPos, d), depth));
	}
	else if (obj.reflT == Refl_t::SPEC)            // Ideal SPECULAR reflection 
	{
		return obj.emi + hitColor.mult(radiance(Ray(hitPos, ray.dir - hitNormal * 2 * hitNormal.dot(ray.dir)), depth));
	}

	Ray reflRay(hitPos, ray.dir - hitNormal * 2 * hitNormal.dot(ray.dir)); // Ideal dielectric REFRACTION 
	bool into = hitNormal.dot(hitNormalNegative) > 0; // Ray from outside going in? 
	double nc = 1;
	double nt = 1.5;
	double nnt = into ? nc / nt : nt / nc;
	double ddn = ray.dir.dot(hitNormalNegative);
	double cos2t;

	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)    // Total internal reflection 
	{
		return obj.emi + hitColor.mult(radiance(reflRay, depth));
	}

	Vec tdir = (ray.dir * nnt - hitNormal * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
	double a = nt - nc, b = nt + nc;
	double R0 = a * a / (b * b);
	double c = 1 - (into ? -ddn : tdir.dot(hitNormal));
	double Re = R0 + (1 - R0) * c * c * c * c * c;
	double Tr = 1 - Re;
	double P = 0.25 + 0.5 * Re;
	double RP = Re / P;
	double TP = Tr / (1 - P);

	return obj.emi + hitColor.mult(depth > 2 ? (erand() < P ?   // Russian roulette 
		radiance(reflRay, depth) * RP : radiance(Ray(hitPos, tdir), depth) * TP) :
		radiance(reflRay, depth) * Re + radiance(Ray(hitPos, tdir), depth) * Tr);
}

int main() {
	int width = 1024, height = 768;
	int samples = Samples <= 4 ? 1 : Samples / 4 ; // # samples 
	Ray camera(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir 
	Vec cx = Vec(width * .5135 / height);
	Vec cy = (cx % camera.dir).norm() * .5135;
	Vec blockRadiance;
	Vec* pixels = new Vec[width * height];

#pragma omp parallel for schedule(dynamic, 1) private(blockRadiance)       // OpenMP 
	for (int y = 0; y < height; y++) // Loop over image rows 
	{
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples * 4, 100. * y / (height - 1.0));
		for (unsigned int x = 0; x < width; x++) // Loop cols 
		{
			for (int subpixelY = 0, i = (height - y - 1) * width + x; subpixelY < 2; subpixelY++) // 2x2 subpixel rows 
			{
				for (int subpixelX = 0; subpixelX < 2; subpixelX++, blockRadiance = Vec()) // 2x2 subpixel cols 
				{
					for (int sample = 0; sample < samples; sample++)
					{
						double r1 = 2 * erand();
						double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);

						double r2 = 2 * erand();
						double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						Vec rayDirection = cx * (((subpixelX + 0.5 + dx) / 2 + x) / width - .5) + cy * (((subpixelY + .5 + dy) / 2 + y) / height - .5) + camera.dir;

						blockRadiance = blockRadiance + radiance(Ray(camera.origin + rayDirection * 140, rayDirection.norm()), 0) * (1.0 / samples);
					} // Camera rays are pushed ^^^^^ forward to start in interior 
					pixels[i] = pixels[i] + Vec(clamp(blockRadiance.x), clamp(blockRadiance.y), clamp(blockRadiance.z)) * .25;
				}
			}
		}
	}

	output_ppm("SmallptVS_Output.ppm", width, height, pixels);

	return (0);
}