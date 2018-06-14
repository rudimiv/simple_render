#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt 
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2 
#include <iostream>

using namespace std;

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm 
	double x, y, z;                  // position, also color (r,g,b) 
	Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; } 

	Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); } 
	Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); } 
	Vec operator*(double b) const { return Vec(x*b,y*b,z*b); } 
	Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }

	Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }

	double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } 
	// cross: 
	Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}

	friend ostream& operator<<(ostream& os, const Vec& v);
}; 

ostream& operator<<(ostream& os, const Vec& v)  
{
    os  << '(' << v.x << ", " << v.y << ", " << v.z << ')';
    return os;  
}  
  

struct Ray {
	Vec start, direction; 
	Ray(Vec o_, Vec d_) : start(o_), direction(d_) {} 
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

		// Solve quadratic equation t^2 * d.direction + 2 * t * (o-p)* d.direction + (o-p) * (o-p)-R^2 = 0 
		//Vec op = p - r.start;
		Vec op = r.start - p;

		double t, eps=1e-4;

		double b = op.dot(r.direction);
		double a = 1.0; // as direction is normed
		double c = op.dot(op) - rad * rad;


		// determinant
		double det = b * b - a * c; 

		if (det < 0) 
			return 0; 
		else 
			det=sqrt(det); 


		return (t = b - det) > eps ? t : ((t=b+det)>eps ? t : 0); 
	}

	Vec normal(const Vec &intersect_point) const {
		return  intersect_point - p;
	}
};

struct Cube : Object {
	double size;
	Vec p, e, c;
	Refl_t refl;
};

struct Plane : Object {
	double p, e, c; // position, emission, color
};




/*void init_scene() {

}*/

/*
	green 008000
	spring green 00FF7F
*/
// массив, работающий с помощтю виртуализации
Object* objects[] = {//Scene: radius, position, emission, color, material 
	new Sphere(Vec(),Vec(.75,.25,.25),DIFF, Vec( 1e5+1,40.8,81.6), 1e5),//Left 
	new Sphere(Vec(),Vec(.25,.25,.75),DIFF, Vec(-1e5+99,40.8,81.6), 1e5),//Rght 
	new Sphere(Vec(),Vec(.75,.75,.75),DIFF, Vec(50,40.8, 1e5), 1e5),//Back 
	new Sphere(Vec(),Vec(), DIFF, Vec(50,40.8,-1e5+170), 1e5),//Frnt 
	new Sphere(Vec(),Vec(.75,.75,.75),DIFF, Vec(50, 1e5, 81.6), 1e5),//Botm 
	new Sphere(Vec(), Vec(.75,.75,.75), DIFF, Vec(50,-1e5+81.6,81.6), 1e5),//Top 
	new Sphere(Vec(),Vec(1,1,1)*.999, DIFF, Vec(27,16.5,47), 16.5),//Mirr 
	new Sphere(Vec(),Vec(1,1,1)*.999, DIFF, Vec(73,16.5,78), 16.5),//Glas 
	new Sphere(Vec(12,12,12),  Vec(), DIFF, Vec(50,681.6-.27,81.6), 600) //Lite */
};

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); } 


// Finds nearest object intersection
inline bool intersect(const Ray &r, double &t, int &id){ 

	// double n = sizeof(objects)/sizeof(Sphere);
	double d;

	double inf=t=1e20;

	double eps=1e-4;
	for(int i=8; i>=0; i--) {
		// std::cout << i << std::endl;
		if(((d = objects[i]->intersect(r)) > 0.5) && d < t) {
			t=d;
			id=i;
		}
	}

	// std::cout << id << " " << t << std::endl;

	return t < inf; 
} 

Vec radiance(const Ray &r, int depth, unsigned short *Xi){ 
	cout.precision(3);
	double t;                               // distance to intersection 
	int id=0;  
	
	// id of intersected object 
	if (!intersect(r, t, id)) return Vec(); // if miss, return black 
	
	Object *obj = objects[id];        // the hit object 
	
	// Intersection point coords
	Vec x = r.start+r.direction*t;

	Vec n = obj->normal(x);
	// std::cout << n;



	Vec nl = n.dot(r.direction) < 0 ? n : n*-1;
	Vec f = obj->c; 
	
	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl 
	
	if (++depth>5) {
		if (erand48(Xi)<p) {
			f=f*(1/p);
		} else {
			// std::cout << "obj->e " << obj->e << std::endl;
			return obj->e; //R.R.
		} 
	}
	
	if (obj->refl == DIFF){                  // Ideal DIFFUSE reflection 
		double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2); 
		Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u; 
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();

		/*if (id == 8) {
			std::cout << obj->e << std::endl;
			Vec res = obj->e + f.mult(radiance(Ray(x,d),depth,Xi));
			std::cout << "DIFF d: " << depth << " id: " << id << " res: "<< res << std::endl;
			return res;
		}*/

		Vec res = obj->e + f.mult(radiance(Ray(x,d),depth,Xi));
		//std::cout << "DIFF d: " << depth << " id: " << id << " res: "<< res << obj->e << obj->c << std::endl;

		return res;

	} else if (obj->refl == SPEC) {           // Ideal SPECULAR reflection 
		Vec res = obj->e + f.mult(radiance(Ray(x,r.direction-n*2*n.dot(r.direction)),depth,Xi)); 
		std::cout << "SPEC d: " << depth << " id: " << id<< " res: "<< res << std::endl;
		return res;
	} else {
		Ray reflRay(x, r.direction-n*2*n.dot(r.direction));     // Ideal dielectric REFRACTION 
		bool into = n.dot(nl)>0;                // Ray from outside going in? 
		double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.direction.dot(nl), cos2t;

		if ((cos2t=1-nnt*nnt*(1-ddn*ddn)) < 0)    // Total internal reflection 
			return obj->e + f.mult(radiance(reflRay,depth,Xi));

		Vec tdir = (r.direction*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
		double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
		double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);

		Vec res = obj->e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette 
		 radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) : 
		 radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr); 

		std::cout << "REFR d: " << depth << " id: " << id << " res: "<< res << std::endl;

		return res;
	}
} 



int main(int argc, char *argv[]){
	if (argc < 3) {
		std::cout << "Please enter rendering parametr and result file name" << std::endl;
		return 0;
	}

	int w = 1024;
	int h=768;
	int samps = atoi(argv[1])/4 ; // # samples

	Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 
	Vec cx=Vec(w*.5135/h), cy=(cx%cam.direction).norm()*.5135, r, *c=new Vec[w * h]; 

	#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 
   	for (int y=0; y<h; y++){                       // Loop over image rows 
    	fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); 
    	for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols 
    	for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows 
         for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols 
           for (int s=0; s<samps; s++){ 
             double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1); 
             double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2); 
             Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) + 
                     cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.direction; 
             Vec rad = radiance(Ray(cam.start+d*140,d.norm()),0,Xi)*(1./samps); 
             //std::cout << y << " " << x << " rad:"<< rad << std::endl;
              r = r + rad;

           } // Camera rays are pushed ^^^^^ forward to start in interior 
           c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25; 
         } 
   } 


	FILE *f = fopen(argv[2], "w");         // Write image to PPM file. 
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 

	for (int i = 0; i < w * h; i++)
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
} 