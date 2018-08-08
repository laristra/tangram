/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef RPG_TOOLS_PRIMITIVES_H_
#define RPG_TOOLS_PRIMITIVES_H_

#include <stdlib.h>
#include <cmath>
#include <math.h>
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"
#include "tangram/support/Matrix.h"

/* Helper functions for generation of basic shapes that can be
   used to specify material distributions */

constexpr double PI = std::acos(-1.0);

inline double pow2(double x) { return x*x; }

/*!
 @brief Matrix for rotating a point around axis centered at the origin 
 by a given angle.

 @param[in] axis Axis of rotation
 @param[in] rotation_angle Angle of rotation
 @return Rotation matrix
*/
Tangram::Matrix rotation_matrix(const Tangram::Vector3& axis,
                                double rotation_angle) {
  Tangram::Matrix rot_matrix(3, 3);
  double cosang = cos(rotation_angle);
  double sinang = sin(rotation_angle);

  rot_matrix[0][0] =  cosang         + axis[0]*axis[0]*(1.0 - cosang);
  rot_matrix[0][1] = -axis[2]*sinang + axis[0]*axis[1]*(1.0 - cosang);
  rot_matrix[0][2] =  axis[1]*sinang + axis[0]*axis[2]*(1.0 - cosang);

  rot_matrix[1][0] =  axis[2]*sinang + axis[0]*axis[1]*(1.0 - cosang);
  rot_matrix[1][1] =  cosang         + axis[1]*axis[1]*(1.0 - cosang);
  rot_matrix[1][2] = -axis[0]*sinang + axis[1]*axis[2]*(1.0 - cosang);

  rot_matrix[2][0] = -axis[1]*sinang + axis[0]*axis[2]*(1.0 - cosang);
  rot_matrix[2][1] =  axis[0]*sinang + axis[1]*axis[2]*(1.0 - cosang);
  rot_matrix[2][2] =  cosang         + axis[2]*axis[2]*(1.0 - cosang);
  
  return rot_matrix;
}

/*!
 @brief Generates a circle with a given center and of a given radius in
 a plane defined by the center and the normal.

 @param[in] center Center of the circle
 @param[in] radius Radius of the circle
 @param[in] normal Normal to the plane
 @param[in] nsamples Number of points representing the circle
 @return Points representing the circle which are ordered ccw
*/
std::vector<Tangram::Point3> circle3d(const Tangram::Point3& center,
                                      double radius,
                                      const Tangram::Vector3& normal,
                                      int nsamples) {
  //By default, we start with the i vector aligned with the x-axis
  Tangram::Vector<3> start_vec(1.0, 0.0, 0.0);
  //Then we project that vector on the circle's plane                             
  start_vec = start_vec - Tangram::dot(start_vec, normal)*normal;
  double vnorm = start_vec.norm();
  //If we projection is too small, we restart with the j vector
  //aligned with the y-axis
  if (vnorm < sqrt(std::numeric_limits<double>::epsilon())) {
    start_vec.axis(1);
    start_vec = start_vec - Tangram::dot(start_vec, normal)*normal;
    vnorm = start_vec.norm();
  }
  //Scale to radius length
  start_vec *= (radius/vnorm);

  std::vector<Tangram::Point3> circle_pts;
  circle_pts.reserve(nsamples);
  //Add the starting point
  circle_pts.push_back(center + start_vec);
  //Get the remaining points by rotation
  for (int isample = 0; isample < nsamples - 1; isample++) {
    double rot_angle = 2*PI*(isample + 1)/nsamples;
    Tangram::Matrix rot_matrix = rotation_matrix(normal, rot_angle);

    circle_pts.push_back(center + rot_matrix*start_vec);
  }

  return circle_pts;
}

/*!
 @brief Generates a cog with a given center and of a given radius in
 a plane defined by the center and the normal.

 @param[in] center Center of the cog
 @param[in] inner_radius Radius of the solid part of the cog
 @param[in] outer_radius Radius bounding the extending teeth
 @param[in] normal Normal to the plane
 @param[in] nteeth Number of teeth the cog has
 @return Points representing the cog which are ordered ccw
*/
std::vector<Tangram::Point3> cog3d(const Tangram::Point3& center,
                                   double inner_radius,
                                   double outer_radius,
                                   const Tangram::Vector3& normal,
                                   int nteeth) {
  assert(nteeth > 2);

  //By default, we start with the i vector aligned with the x-axis
  Tangram::Vector<3> outer_start_vec(1.0, 0.0, 0.0);
  //Then we project that vector on the circle's plane                             
  outer_start_vec = outer_start_vec - Tangram::dot(outer_start_vec, normal)*normal;
  double vnorm = outer_start_vec.norm();
  //If we projection is too small, we restart with the j vector
  //aligned with the y-axis
  if (pow2(vnorm) < std::numeric_limits<double>::epsilon()) {
    outer_start_vec.axis(1);
    outer_start_vec = outer_start_vec - Tangram::dot(outer_start_vec, normal)*normal;
    vnorm = outer_start_vec.norm();
  }
  Tangram::Vector<3> inner_start_vec = outer_start_vec;
  //Scale to the espective radius length
  outer_start_vec *= (outer_radius/vnorm);
  inner_start_vec *= (inner_radius/vnorm);

  std::vector<Tangram::Point3> cog_pts;
  cog_pts.reserve(4*nteeth);

  //We want chords on the inner and the outer circles to be of the same length,
  //and the segments between circles to be of the same length
  double tooth_angle = 2*PI/nteeth;
  double inner_angle = 2*std::atan(
    sin(tooth_angle/4.0)/(inner_radius/outer_radius + cos(tooth_angle/4.0)) );
  double dtooth_ang[4] = { 0.0, tooth_angle/2.0 - inner_angle,
                           tooth_angle/4.0, inner_angle };

  //Get the points by rotation
  for (int itooth = 0; itooth < nteeth; itooth++) {
    double rot_angle = 2*PI*itooth/nteeth;
    for (int its = 0; its < 4; its++) {
      rot_angle += dtooth_ang[its];
      Tangram::Matrix rot_matrix = rotation_matrix(normal, rot_angle);

      if (its < 2)
        cog_pts.push_back(center + rot_matrix*outer_start_vec);
      else
        cog_pts.push_back(center + rot_matrix*inner_start_vec);
    }
  }

  return cog_pts;
}

/*!
 @brief Generates a prism given one of its bases. The prism is
 extruded in the direction opposite to the base's normal.

 @param[in] base_pts Points defining the base of the prism
 @param[in] height Height of the extruded prism
 @param[in] base_scaling Scaling factor for the opposite face; 
 if set to 1 the prism will be a right prism
 @param[in] base_normal Normal to the given base; if not
 specified, it will be computed
 @return MatPoly object corresponding to the resulting prism
*/
Tangram::MatPoly<3> prism(const std::vector<Tangram::Point3>& base_pts,
                          double height,
                          double base_scaling,
                          const Tangram::Vector3& base_normal = 
                            Tangram::Vector3(0.0, 0.0, 0.0)) {
  assert(base_scaling > std::numeric_limits<double>::epsilon());

  int nbase_vrts = (int) base_pts.size();
  std::vector< std::vector<int> > ifaces_vrts(nbase_vrts + 2);
  ifaces_vrts[0].resize(nbase_vrts);
  std::iota(ifaces_vrts[0].begin(), ifaces_vrts[0].end(), 0);

  std::vector<Tangram::Point3> prism_pts;
  prism_pts.reserve(2*nbase_vrts);
  prism_pts.insert(prism_pts.begin(), base_pts.begin(), base_pts.end());

  Tangram::Vector3 base_shift;
  if (base_normal.is_zero())
   base_shift = Tangram::polygon3d_normal(prism_pts, ifaces_vrts[0]);
  else
   base_shift = base_normal;

  base_shift *= -height;     

  for(int ivrt = 0; ivrt < nbase_vrts; ivrt++)
    prism_pts.push_back(Tangram::Point3(
      prism_pts[(nbase_vrts - ivrt)%nbase_vrts] + base_shift));

  for(int ivrt = 0; ivrt < nbase_vrts; ivrt++) {
    int iface = ivrt + 1;
    ifaces_vrts[iface] = { (ivrt + 1)%nbase_vrts, ivrt,
      nbase_vrts + (nbase_vrts - ivrt)%nbase_vrts,
      2*nbase_vrts - ivrt - 1 };
  }

  ifaces_vrts[nbase_vrts + 1].resize(nbase_vrts);
  std::iota(ifaces_vrts[nbase_vrts + 1].begin(), 
            ifaces_vrts[nbase_vrts + 1].end(), nbase_vrts);

  if (base_scaling != 1.0) {    
    std::vector<double> new_base_moments;
    Tangram::polygon3d_moments(prism_pts, ifaces_vrts[nbase_vrts + 1], 
                               new_base_moments);

    Tangram::Point3 new_base_cen;
    for (int idim = 0; idim < 3; idim++)
      new_base_cen[idim] = new_base_moments[idim + 1]/new_base_moments[0];

    for(int ivrt = 0; ivrt < nbase_vrts; ivrt++) {
      Tangram::Vector3 cen2vrt = prism_pts[nbase_vrts + ivrt] - new_base_cen;
      cen2vrt *= base_scaling;
      prism_pts[nbase_vrts + ivrt] = new_base_cen + cen2vrt;
    }
  }

  Tangram::MatPoly<3> mat_poly;
  mat_poly.initialize(prism_pts, ifaces_vrts);

  return mat_poly;  
}

/*!
 @brief Generates a prism with non-parallel bases. The prism
 can become self-intersecting if overly twisted.

 @param[in] base_pts Points defining the base of the prism
 @param[in] opp_base_normal Normal to the opposite base
 @param[in] opp_base_centroid Centroid of the opposite base
 @param[in] base_scaling Scaling factor for the opposite face; 
 if set to 1 the bases will be the same wrt to their respective planes
 @param[in] base_normal Normal to the given base; if not
 specified, it will be computed
 @return MatPoly object corresponding to the resulting prism
*/
Tangram::MatPoly<3> skewed_prism(const std::vector<Tangram::Point3>& base_pts,
                                 const Tangram::Vector3& opp_base_normal,
                                 const Tangram::Point3& opp_base_centroid,
                                 double base_scaling,
                                 const Tangram::Vector3& base_normal = 
                                    Tangram::Vector3(0.0, 0.0, 0.0)) {
  assert(base_scaling > sqrt(std::numeric_limits<double>::epsilon()));

  int nbase_vrts = (int) base_pts.size();
  std::vector< std::vector<int> > ifaces_vrts(nbase_vrts + 2);
  ifaces_vrts[0].resize(nbase_vrts);
  std::iota(ifaces_vrts[0].begin(), ifaces_vrts[0].end(), 0);

  std::vector<Tangram::Point3> prism_pts;
  prism_pts.reserve(2*nbase_vrts);
  prism_pts.insert(prism_pts.begin(), base_pts.begin(), base_pts.end());

  std::vector<double> base_moments;
  Tangram::polygon3d_moments(prism_pts, ifaces_vrts[0], base_moments);

  Tangram::Point3 base_centroid;
  for (int idim = 0; idim < 3; idim++)
    base_centroid[idim] = base_moments[idim + 1]/base_moments[0];

  Tangram::Vector3 base_start_vec = prism_pts[0] - base_centroid;
  double vnorm = base_start_vec.norm();
  assert(vnorm > sqrt(std::numeric_limits<double>::epsilon()));
  base_start_vec /= vnorm;

  //We find the starting point on the opposite plane by intersecing it
  //with the line parallel to the line conneting centroids
  Tangram::Vector3 line_dvec = opp_base_centroid - base_centroid;
  vnorm = line_dvec.norm();
  assert(vnorm > sqrt(std::numeric_limits<double>::epsilon()));
  line_dvec /= vnorm;

  double denom = Tangram::dot(line_dvec, opp_base_normal);
  assert(denom > sqrt(std::numeric_limits<double>::epsilon()));

  Tangram::Vector3 opp_base_start_vec = prism_pts[0] - opp_base_centroid +
    (Tangram::dot(opp_base_centroid - prism_pts[0], opp_base_normal)/denom)*line_dvec;

  vnorm = opp_base_start_vec.norm();
  assert(vnorm > sqrt(std::numeric_limits<double>::epsilon()));
  opp_base_start_vec /= vnorm; 

  //Base normal to check the rotation angle sign
  Tangram::Vector3 base_normal_vec;
  if (base_normal.is_zero())
    base_normal_vec = Tangram::polygon3d_normal(prism_pts, ifaces_vrts[0]);
  else 
    base_normal_vec = base_normal;

  //We find opposite base vertices by rotating them by the same angle
  for(int ivrt = 0; ivrt < nbase_vrts; ivrt++) {
    int icur_vrt = (nbase_vrts - ivrt)%nbase_vrts;

    Tangram::Vector3 base_cen2vrt_vec = prism_pts[icur_vrt] - base_centroid;
    double base_cen2vrt_dst = base_cen2vrt_vec.norm();
    base_cen2vrt_vec /= base_cen2vrt_dst;

    //acos gives the value from 0 to PI, need sin to determine the angle sign
    int sign_rot_angle = Tangram::dot(Tangram::cross(base_start_vec, base_cen2vrt_vec), 
                                      base_normal_vec) > 0.0 ? 1 : -1;
    
    double cos_ang = Tangram::dot(base_cen2vrt_vec, base_start_vec);
    //Due to numerical errors, the value can go out of the cos range
    if (cos_ang > 1.0 - std::numeric_limits<double>::epsilon()) cos_ang = 1.0;
    if (cos_ang < -1.0 + std::numeric_limits<double>::epsilon()) cos_ang = -1.0;

    double rot_angle = sign_rot_angle*std::acos(cos_ang);

    //Normals are outward, so we flip the direction of rotation
    Tangram::Matrix rot_matrix = rotation_matrix(opp_base_normal, -rot_angle);

    //Scale and get vertex coordinate
    prism_pts.push_back(opp_base_centroid + 
      base_scaling*base_cen2vrt_dst*(rot_matrix*opp_base_start_vec));
  }

  for(int ivrt = 0; ivrt < nbase_vrts; ivrt++) {
    int iface = ivrt + 1;
    ifaces_vrts[iface] = { (ivrt + 1)%nbase_vrts, ivrt,
      nbase_vrts + (nbase_vrts - ivrt)%nbase_vrts,
      2*nbase_vrts - ivrt - 1 };
  }

  ifaces_vrts[nbase_vrts + 1].resize(nbase_vrts);
  std::iota(ifaces_vrts[nbase_vrts + 1].begin(), 
            ifaces_vrts[nbase_vrts + 1].end(), nbase_vrts);

  std::vector<double> prj_poly_moments;
  Tangram::polygon3d_moments(prism_pts, ifaces_vrts[nbase_vrts + 1], 
                              prj_poly_moments);

  Tangram::MatPoly<3> mat_poly;
  mat_poly.initialize(prism_pts, ifaces_vrts);

  return mat_poly;  
}

/*!
 @brief Generates a faceted sphere.

 @param[in] center Center of the sphere
 @param[in] radius Radius of the sphere
 @param[in] nquadrant_samples For each of the polar angles, the 
 number of values sampled per a circle's quadrant
 @return MatPoly object corresponding to the resulting sphere;
 the number of faces is 8*pow2(nquadrant_samples)
*/
Tangram::MatPoly<3> sphere(const Tangram::Point3& center,
                           double radius,
                           int nquadrant_samples) {
  std::vector<Tangram::Point3> sphere_pts;
  int ncircle_samples = 4*nquadrant_samples;
  sphere_pts.reserve((2*nquadrant_samples - 1)*ncircle_samples + 2);

  //Top point
  sphere_pts.push_back(Tangram::Point3(center[0], center[1], center[2] + radius));
  //Level circles
  for (int iin = 0; iin < 2*nquadrant_samples - 1; iin++) {
    double inclination = (iin + 1)*PI/(2*nquadrant_samples);
    for (int iaz = 0; iaz < ncircle_samples; iaz++) {
      double azimuth = iaz*PI/(2*nquadrant_samples);
      sphere_pts.push_back(Tangram::Point3(
        center[0] + radius*sin(inclination)*cos(azimuth),
        center[1] + radius*sin(inclination)*sin(azimuth),
        center[2] + radius*cos(inclination)));
    }
  }
  //Bottom point
  sphere_pts.push_back(Tangram::Point3(center[0], center[1], center[2] - radius));

  std::vector< std::vector<int> > ifaces_vrts(8*pow2(nquadrant_samples));
  //Top triangular faces
  for (int iaz = 0; iaz < ncircle_samples; iaz++)
    ifaces_vrts[iaz] = {0, 1 + iaz, 1 + (iaz + 1)%ncircle_samples};

  //Quad faces
  for (int iin = 0; iin < 2*(nquadrant_samples - 1); iin++) {
    int offset = (iin + 1)*ncircle_samples;
    for (int iaz = 0; iaz < ncircle_samples; iaz++)
      ifaces_vrts[offset + iaz] = {
        offset - ncircle_samples + 1 + (iaz + 1)%ncircle_samples,
        offset - ncircle_samples + 1 + iaz, 
        offset + 1 + iaz,
        offset + 1 + (iaz + 1)%ncircle_samples };
  }

  //Bottom triangular faces
  int offset = (2*nquadrant_samples - 1)*ncircle_samples;
  for (int iaz = 0; iaz < ncircle_samples; iaz++)
    ifaces_vrts[offset + iaz] = {offset + 1, 
      offset - ncircle_samples + 1 + (iaz + 1)%ncircle_samples, 
      offset - ncircle_samples + 1 + iaz };

  Tangram::MatPoly<3> mat_poly;
  mat_poly.initialize(sphere_pts, ifaces_vrts);

  return mat_poly;      
}

#endif