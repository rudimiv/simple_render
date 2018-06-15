#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2 
#include <iostream>
#include <vector>

using namespace std;

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm 
	double x, y, z;                  // position, also color (r,g,b) 
	Vec(double x_=0, double y_=0, double z_=0){ x = x_; y = y_; z = z_; } 
	Vec operator+(const Vec &b) const { return Vec(x + b.x,y + b.y,z + b.z); } 
	Vec operator-(const Vec &b) const { return Vec(x - b.x,y - b.y,z - b.z); } 
	Vec operator*(double b) const { return Vec(x * b,y * b,z * b); } 
	Vec mult(const Vec &b) const { return Vec(x * b.x,y * b.y,z * b.z); } 
	Vec& norm(){ return *this = *this * (1/sqrt(x * x + y * y + z * z)); } 
	double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross
	double length() const {return sqrt(x * x + y * y + z * z);}
	Vec operator%(const Vec&b) const {return Vec(y * b.z - z * b.y,z * b.x - x * b.z,x * b.y - y * b.x);} 

	friend ostream& operator<<(ostream& os, const Vec& v);
}; 

ostream& operator<<(ostream& os, const Vec& v)  
{
		os  << '(' << v.x << ", " << v.y << ", " << v.z << ')';
		return os;  
} 

struct Ray {
	Vec o, d; 
	Ray(Vec o_, Vec d_) : o(o_), d(d_) {} 
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance() 




struct Object {
	// returns distance, 0 if nohit
	virtual double intersect(const Ray &r) const = 0;
	virtual Vec normal(const Vec &intersect_point) const = 0;

	Object(Vec e_, Vec c_, Refl_t refl_): e(e_), c(c_), refl(refl_) {}
	
	Vec e, c;	// emission, color
	Refl_t refl;	// reflection type (DIFFuse, SPECular, REFRactive)
};

struct Sphere : public Object { 
	double rad;	// radius 
	Vec p;	// position 

	Sphere(Vec e_, Vec c_, Refl_t refl_, Vec p_, double rad_):  Object(e_, c_, refl_), rad(rad_), p(p_) {} 

	double intersect(const Ray &r) const {

		// Solve quadratic equation t^2 * d.direction + 2 * t * (o - p)* d.direction + (o - p) * (o - p)-R^2 = 0 
		Vec op = p - r.o;
		// Vec op = r.start - p;

		double t, eps = 1e-4;

		double b = op.dot(r.d);
		double a = 1.0; // as direction is normed
		double c = op.dot(op) - rad * rad;


		// determinant
		double det = b * b - a * c; 

		if (det < 0) 
			return 0; 
		else 
			det = sqrt(det); 


		return (t = b - det) > eps ? t : ((t = b + det)>eps ? t : 0); 
	}

	Vec normal(const Vec &intersect_point) const {
		return  (intersect_point - p).norm();
	}
};

struct Rectangle : Object {
	Vec p; // left corner
	Vec vector_a, vector_b; // borders vectors
	double len_a, len_b;
	Vec n;

	Rectangle(Vec e_, Vec c_, Refl_t refl_, Vec A, Vec B, Vec C): Object(e_, c_, refl_),
																p(C), 
																vector_a(A - C),
																vector_b(B - C) { 
		n = (vector_a % vector_b).norm();
		len_a = vector_a.length();
		len_b = vector_b.length();
	}

	double intersect(const Ray &r) const {
		// Check intersection with plane
		double nv = n.dot(r.d);
		double eps = 1e-4;

		// std::cout << "*****************\n";
		// std::cout << "n: " << n << std::endl;
		// std::cout << "nv: " << nv << std::endl;

		if (fabs(nv) < eps) {
			return 0;
		}

		// std::cout << "vector_a: " << vector_a << std::endl;
		// std::cout << "vector_b: " << vector_b << std::endl;
		// std::cout << "dir_vector: " <<  "o: " << r.o << " d: " << r.d << std::endl;
		// std::cout << "p: " << p << std::endl;

		Vec po = p - r.o;
		
		// std::cout << "po: " << po << std::endl;
		double k = n.dot(po) / nv;
		
		// Intersection point
		Vec I = r.o + r.d * k;

		// std::cout << "k: " << k << std::endl;
		// std::cout << "I: " << I << std::endl;
		// Does Intersection point belong to rectangle?
		Vec ip = I - p;

		// triangle square
		// Vec res = vector_a % vector_b;
		double square = (vector_a % vector_b).length() / 2;
		double square_1 = (ip % vector_a).length() / 2;
		double square_2 = (ip % vector_b).length() / 2;
		double square_3 = ((ip - vector_b) % (vector_b - vector_a)).length() / 2;

		//std::cout << "ip: " << ip << std::endl;

		/*if (I.x >= 20 && I.x <= 50 && I.y <= 20 && I.y >= 0){
			std::cout << "I: " << I << std::endl;
			std::cout << square - square_1 - square_2 - square_3 << std::endl;
		}*/
		if (fabs(square - square_1 - square_2 - square_3 ) > eps) {
			return 0;
		}

		// std::cout << "K_res: " << k << std::endl;

		// std::cout << "I: " << I << std::endl;
		return k;
	}

	Vec normal(const Vec &intersect_point) const {
		return n;
	}
};


std::vector < Object*> objects;

void init_scene() {
	uint pyramid_shift_y = 60 + 20;
	uint pyramid_shift_x = 10;

	Vec green = Vec(0,.501,0);
	Vec spring_green = Vec(0,1,.49);
	Vec red = Vec(1, 0, 0);
	Vec gray = Vec(.75,.75,.75);
	Vec blue = Vec(0,0,1);
	Vec white = Vec(1,1,1);

	//Scene: emission, color, material, position, radius
	objects.push_back(new Sphere(Vec(), green ,DIFF, Vec( 1e5 + 1,40.8,81.6), 1e5)); //Left 
	objects.push_back(new Sphere(Vec(), red,DIFF, Vec(-1e5 + 99,40.8,81.6), 1e5)); //Rght 
	objects.push_back(new Sphere(Vec(), blue,DIFF, Vec(50,40.8, 1e5), 1e5)); //Back 
	objects.push_back(new Sphere(Vec(),Vec(), DIFF, Vec(50,40.8,-1e5 + 170), 1e5)); //Frnt 
	objects.push_back(new Sphere(Vec(), white ,DIFF, Vec(50, 1e5, 81.6), 1e5)); //Botm 
	objects.push_back(new Sphere(Vec(), Vec(), DIFF, Vec(50,-1e5 + 81.6,81.6), 1e5)); //Top 
	objects.push_back(new Sphere(Vec(),blue*.999, REFR, Vec(27,16.5,47), 16.5)); //Mirr 
	objects.push_back(new Sphere(Vec(),Vec(1,1,1)*.999, SPEC, Vec(73,16.5,78), 16.5)); //Glas 
	objects.push_back(new Sphere(Vec(12,12,12),  Vec(), DIFF, Vec(50,681.6-.27,81.6), 600)); //Lite
	

	Vec A_point = Vec(20 + pyramid_shift_x, 0, pyramid_shift_y);
	Vec B_point = Vec(50 + pyramid_shift_x, 0, pyramid_shift_y);
	Vec C_point = Vec(35 + pyramid_shift_x,0,pyramid_shift_y + 26);
	Vec D_point = Vec(35 + pyramid_shift_x, 30 * 0.81, pyramid_shift_y + 8.66);
	// Pyramids
	objects.push_back(new Rectangle(Vec(), Vec(1,1,1)*.999, SPEC, A_point, B_point, D_point));
	objects.push_back(new Rectangle(Vec(), Vec(1,1,1)*.999, SPEC, B_point, C_point, D_point));
	objects.push_back(new Rectangle(Vec(), Vec(1,1,1)*.999, SPEC, A_point, C_point, D_point));
	
}



inline double clamp(double x){ return x < 0 ? 0 : x > 1 ? 1 : x; } 
inline int toInt(double x){ return int(pow(clamp(x),1/2.2) * 255 + .5); } 


// Finds nearest object intersection
inline bool intersect(const Ray &r, double &t, int &id){ 

	// double n = sizeof(objects)/sizeof(Sphere);
	double d;

	double inf = t = 1e20;

	double eps = 1e-4;
	for(int i = objects.size() - 1; i >= 0; i--) {
		// std::cout << i << std::endl;
		if((d = objects[i]->intersect(r)) && d < t && d >= 0) {
			t = d;
			id = i;
		}
	}

	// std::cout << id << " " << t << std::endl;
	/*if (id == objects.size() - 1) {
		std::cout << "yes " << t << std::endl;
	} else {
		std::cout << "no" << t << std::endl;
	}*/
	return t < inf; 
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi){ 
	double t;                               // distance to intersection 

	int id = 0;                               // id of intersected object 

	if (!intersect(r, t, id)) {
		return Vec(); // if miss, return black 
	}

	const Object *obj = objects[id];        // the hit object 
	
	

	// Intersection point coords
	Vec x = r.o + r.d * t;

	// Normal to object
	Vec n = obj->normal(x);
	//Vec n=(x - obj->p).norm();

	Vec nl = n.dot(r.d) < 0 ? n : n * -1;
	Vec f = obj->c; 


	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl 

	if (++depth > 5) {
		if (erand48(Xi) < p) {
			f = f * (1/p); 
		} else {
			return obj->e; //R.R.
		}
	}

	if (obj->refl == DIFF){
		// Ideal DIFFUSE reflection 
		double r1 = 2 * M_PI * erand48(Xi);
		double r2 = erand48(Xi);
		double r2s = sqrt(r2); 

		Vec w = nl;
		Vec u = ((fabs(w.x)>.1 ? Vec(0,1, 0) : Vec(1, 0, 0)) % w).norm() ; 
		Vec v = w % u;

		// сos распределение
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
		Vec res = obj->e + f.mult(radiance(Ray(x,d),depth,Xi));
		//std::cout << "DIFF d: " << depth << " id: " << id << " res: "<< res << obj.e << obj.c << std::endl;
		return res;
	} else if (obj->refl == SPEC) {        
		// Ideal SPECULAR reflection 
		return obj->e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi)); 
	} else {
		// Ideal dielectric REFRACTION 

		Ray reflRay(x, r.d - n * 2 * n.dot(r.d));
		
		bool into = n.dot(nl) > 0;                // Ray from outside going in? 
		
		double nc = 1, nt = 1.5, nnt = into?nc/nt : nt/nc, ddn = r.d.dot(nl), cos2t; 
		
		if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn))<0)    // Total internal reflection 
			return obj->e + f.mult(radiance(reflRay,depth,Xi)); 
		
		Vec tdir = (r.d * nnt - n*((into?1:-1)*(ddn * nnt + sqrt(cos2t)))).norm(); 
		
		double a = nt - nc, b = nt + nc, R0 = a * a/(b * b), c = 1-(into?-ddn : tdir.dot(n)); 
		double Re=R0+(1 - R0)*c * c * c * c * c,Tr = 1 - Re,P=.25+.5 * Re,RP = Re/P,TP = Tr/(1 - P); 
		
		return obj->e + f.mult(depth > 2 ? (erand48(Xi)<P ?   // Russian roulette 
			radiance(reflRay,depth,Xi)*RP : radiance(Ray(x,tdir),depth,Xi)*TP) : 
			radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 
	}
} 


int main(int argc, char *argv[]){ 
	if (argc < 3) {
		std::cout << "Please enter rendering parametr and result file name" << std::endl;
		return 0;
	}

	init_scene();

	//int w = 300;
	//int h = 210;

	int w = 1024;
	int h = 768;
	int samps = atoi(argv[1])/4 ; // # samples
	std::cout << samps << std::endl;

	Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 
	Vec cx = Vec(w * .5135/h);
	Vec cy = (cx % cam.d).norm() * .5135;
	Vec r;
	Vec *c = new Vec[w * h]; 
	#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 
	for (int y = 0; y < h; y++){                       // Loop over image rows 
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps * 4,100.*y/(h - 1)); 
		for (unsigned short x = 0, Xi[3]={0,0,y * y * y}; x < w; x++)   // Loop cols 
			for (int sy = 0, i=(h - y - 1)*w + x; sy < 2; sy++)     // 2x2 subpixel rows 
			for (int sx = 0; sx < 2; sx++, r = Vec()){        // 2x2 subpixel cols 
				for (int s = 0; s < samps; s++){ 
					double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1)-1: 1 - sqrt(2 - r1); 
					double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2)-1: 1 - sqrt(2 - r2); 
					Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) + 
								 cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d; 
					Vec rad = radiance(Ray(cam.o + d * 140,d.norm()),0,Xi)*(1./samps); 
					r = r + rad;

				} // Camera rays are pushed ^^^^^ forward to start in interior 
				c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25; 
			} 
	} 

	FILE *f = fopen(argv[2], "w");         // Write image to PPM file. 
	fprintf(f, "P3\n %d %d\n %d\n", w, h, 255); 
	for (int i = 0; i < w * h; i++) 
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
 } 