 #ifndef ARRAY_COMPUTATION_H
 #define ARRAY_COMPUTATION_H

 #include <cmath>

 inline void Substract(double result[2], const double a[2], const double b[2]) {
     result[0] = a[0] - b[0];
     result[1] = a[1] - b[1];
 }

 inline double Norm(const double v[2]) {
     return sqrt(v[0] * v[0] + v[1] * v[1]);
 }

 void periodic_boundary(double mir_pos[],double org_pos[], int adj[]){
    for(int i=0;i<2;i++){
        mir_pos[i]=org_pos[i]+adj[i];
	}
}

 inline double wrap_box(double x) {
     x = fmod(x, 1.0);
     if (x < 0.0) x += 1.0;
     return x;
 }

 #endif
