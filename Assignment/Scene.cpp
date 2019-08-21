#include "Scene.h"

#include "Colour.h"
#include "ImageDisplay.h"
#include "utility.h"

Scene::Scene() : backgroundColour(0,0,0), ambientLight(0,0,0), maxRayDepth(3), renderWidth(800), renderHeight(600), filename("render.png"), camera_(), objects_(), lights_() {

}

Scene::~Scene() {

}

void Scene::render() const {
	ImageDisplay display("Render", renderWidth, renderHeight);
	



	const double w = double(renderWidth);
	const double h = double(renderHeight);

	for (unsigned int v = 0; v < renderHeight; ++v) {
		for (unsigned int u = 0; u < renderWidth; ++u) {
			double cu = -1 + (u + 0.5)*(2.0 / w);
			double cv = -h/w + (v + 0.5)*(2.0 / w);
			Ray ray = camera_->castRay(cu, cv);
			display.set(u, v, computeColour(ray, maxRayDepth));
		}
		display.refresh();
	}

	display.save(filename);
	display.pause(5);
}

RayIntersection Scene::intersect(const Ray& ray) const {
	RayIntersection firstHit;
	firstHit.distance = infinity;	
	for (const auto & obj : objects_) {
		for (const auto & hit : obj->intersect(ray)) {
			if (hit.distance > epsilon && hit.distance < firstHit.distance) {
				firstHit = hit;
			}
		}
	}	return firstHit;
}

Colour Scene::computeColour(const Ray& ray, unsigned int rayDepth) const {
    RayIntersection hitPoint = intersect(ray);
    if (hitPoint.distance == infinity) {
        return backgroundColour;
    }

    Colour hitColour(0, 0, 0);
		
    // Code to do better lighting, shadows, and reflections goes here.
    for (const auto & light: lights_) {
        // Compute the influence of this light on the appearance of the hit object.
        if (light->getDistanceToLight(hitPoint.point) < 0) {
            // An ambient light, ignore shadows and add appropriate colour
            hitColour += light->getIlluminationAt(hitPoint.point) * hitPoint.material.ambientColour;
        } else {

            /*************************************************
             * TODO - ADD DIFFUSE AND SPECULAR LIGHTING HERE *
             *      - SHADOW COMPUTATIONS GO HERE ALSO       *
             *************************************************/
            //Getting the illumination with diffuse and specular colour
            Colour i = light -> getIlluminationAt(hitPoint.point);
            Colour kd = hitPoint.material.diffuseColour;
            Colour ks = hitPoint.material.specularColour;

            //Surface normal
            Vector n = hitPoint.normal;
            //vector from surface point to light source
            Vector l = -(light -> getLightDirection(hitPoint.point));
            n = n/n.norm();
            l = l/l.norm();
            //reflection of l about n
            Vector r = 2 * n.dot(l) * n - l;           
            r = r/r.norm();
            //surface point to view direction
            Vector v = -ray.direction;
            v = v/v.norm();

            //specular exponent
            double specEx = hitPoint.material.specularExponent;
            //Dot product results, making sure its positive or 0
            double nDotL = std::max<double>(0, n.dot(l));
            double rDotV = std::max<double>(0, r.dot(v));

            //Shadow rays
            Ray shadow;
            shadow.direction = l;
            shadow.point = hitPoint.point;
            RayIntersection shadowPoint = intersect(shadow);

            //Adding the light/colour to the object if it is not in a shadow
            if(shadowPoint.distance > light->getDistanceToLight(hitPoint.point)){
                hitColour += i * ((kd * nDotL) + (ks * pow(rDotV, specEx)));
            }

                    

                    
        }
    }

    /**********************************
     * TODO - ADD MIRROR EFFECTS HERE *
     **********************************/

    //Implementing mirror reflections recursively 
    if(rayDepth > 0){
        Ray mirrorRay;

        Vector n = hitPoint.normal;
        Vector v = -ray.direction;
        n = n/n.norm();
        v = v/v.norm();

        //Similar to shadow ray except the direction is the angle of v about n 
        mirrorRay.direction = 2 * n.dot(v) * n - v;
        mirrorRay.point = hitPoint.point;

        hitColour += hitPoint.material.mirrorColour * computeColour(mirrorRay, rayDepth - 1);
    }


    hitColour.clip();
    return hitColour;
}

bool Scene::hasCamera() const {
    return bool(camera_);
}
