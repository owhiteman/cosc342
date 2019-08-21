#include "Plane.h"

#include "utility.h"

Plane::Plane() : Object() {

}

Plane::Plane(const Plane& plane) : Object(plane) {

}

Plane::~Plane() {

}

const Plane& Plane::operator=(const Plane& plane) {
    if (this != &plane) {
        Object::operator=(plane);
    }
    return *this;
}

std::vector<RayIntersection> Plane::intersect(const Ray& ray) const {

    //Transforming the ray
    std::vector<RayIntersection> result;
    Ray inverseRay = transform.applyInverse(ray);

    //intersection
    double z0 = inverseRay.point(2);
    double dz = inverseRay.direction(2);

    double t = -z0/dz;

    RayIntersection hit;

    //0 division check
    if(std::abs(dz) < epsilon){
        return result;
    }
        
    if(t > 0){    
        hit.point = inverseRay.point + t*inverseRay.direction;

        //checking the x y value range, then computing the details of the hit point
        if((hit.point(0) > -1 && hit.point(0) < 1) && (hit.point(1) > - 1 && hit.point(1) < 1)){
            hit.material = material;
                
            hit.point = transform.apply(Point(hit.point));
                
            hit.normal = Normal(0, 0, 1);
            hit.normal = transform.apply(Normal(hit.normal));

            hit.distance =  (hit.point - ray.point).norm();
            result.push_back(hit);
        }
    }


    return result;
}
