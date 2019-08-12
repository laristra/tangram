/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef MATDATA_ROTOR_3D_H_
#define MATDATA_ROTOR_3D_H_

#include <stdlib.h>
#include "tangram/support/tangram.h"
#include "tangram/utility/rpgtools/cuts.h"
#include "tangram/utility/rpgtools/primitives.h"

#include "tangram/driver/CellMatPoly.h"
#include "tangram/driver/write_to_gmv.h"

double find_matching_radius(const double tooth_len,
                            const double tooth_depth,
                            const int nteeth) {
  double l = 0.5*tooth_len, h = 0.5*tooth_depth;
  double t = cos(PI/(2*nteeth));

  double rad = sqrt( pow2(h) + 
    ( pow2(l)*(t + 1) + l*sqrt( pow2(l)*pow2(t + 1) + 4*pow2(h)*(1 - pow2(t)) ) ) / 
    ( 1 - pow2(t) ) );

  return rad;
}

/*!
 @brief For a given position and orientation cuts out a planetary gears system from the
 collection of r3d_poly's. 

 @param[in] polys_data Data on subdomain containing the planetaty gears.
 @param[in] normal Defines the orientation of planetary gears.
 @param[in] frame_center Defines the positon of planetary gears.
 @param[in] frame_depth Defines the depth/thickness of the planetary gears system.
 @param[in] sun_shaft_rad Radius of the main shaft where it is attached to the sun gear
 @param[in] sun_in_rad Radius of solid part of the sun gear
 @param[in] sun_out_rad Radius of the sun gear's teeth circumference
 @param[in] carrier_clearance Distance between carriers and frame/gears
 @param[in] carrier_depth Depth/thickness of carriers
 @param[in] frame_rad Radius of the planetary gears frame
 @param[in] nplanets Number of planetary gears
 @param[in] sun_nteeth Number of teeth for the sun gear
 @param[in] planet_nteeth Number of teeth for planetary gears
 @param[in] shaft_nsides Number of sides for in/out shafts
 @param[in] gasket_r_scaling Relative radius of gaskets with respect to shaft
 @param[in] vol_tol Volume tolerance
 @param[in] dst_tol Distance tolerance
 @param[out] res_poly_sets_data Data on subdomains corresponding 
 to different components of the planetary gears system
*/
Tangram::MatPoly<3>
planetary_gear(const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
               const Tangram::Vector3& normal,
               const Tangram::Point3& frame_center,
               const double frame_depth,
               const double sun_shaft_rad,
               const double sun_in_rad,
               const double sun_out_rad,
               const double carrier_clearance,
               const double carrier_depth,
               const double frame_rad,
               const int nplanets,
               const int sun_nteeth,
               const int planet_nteeth,
               const int shaft_nsides,
               const double gasket_r_scaling,
               const double vol_tol,
               const double dst_tol,
               std::vector< std::vector< std::shared_ptr<RefPolyData_t> > >& 
                res_poly_sets_data) {
  assert((sun_nteeth%nplanets == 0) && (planet_nteeth%2 == 0));

  res_poly_sets_data.clear();

  double gears_depth = frame_depth - 4*carrier_clearance - 2*carrier_depth;
  Tangram::Point3 gears_center = frame_center + 
    (-2*carrier_clearance - carrier_depth)*normal;
  Tangram::Point3 carrier_centers[2] = {
    frame_center + (-carrier_clearance)*normal, 
    frame_center + (-frame_depth + carrier_clearance + carrier_depth)*normal };

  std::vector<Tangram::Point3> sun_pts = 
    cog3d(gears_center, sun_in_rad, sun_out_rad, normal, sun_nteeth, dst_tol);    

  double tooth_len = (sun_pts[1] - sun_pts[0]).norm();
  double tooth_depth = sun_out_rad - sun_in_rad;

  double planet_rad = find_matching_radius(tooth_len, tooth_depth, planet_nteeth);

  double min_gap = 0.0;
  int ring_nteeth = sun_nteeth;
  double ring_out_rad = sun_out_rad;

  while (ring_out_rad < sun_in_rad + 2*(planet_rad + min_gap) + tooth_depth) {
    ring_nteeth += nplanets;
    ring_out_rad = find_matching_radius(tooth_len, tooth_depth, ring_nteeth) + 
                   0.5*tooth_depth;
  }

  double gap = 0.5*(ring_out_rad - sun_in_rad - tooth_depth) - planet_rad;

  std::vector<Tangram::Point3> sun_shaft_pts = 
    circle3d(frame_center, sun_shaft_rad, normal, shaft_nsides, dst_tol);
  Tangram::MatPoly<3> sun_shaft_poly = prism(sun_shaft_pts, 
    2*carrier_clearance + carrier_depth + gears_depth, 1.0, dst_tol, normal);

  Tangram::MatPoly<3> sun_poly = prism(sun_pts, gears_depth, 1.0, dst_tol, normal);

  std::vector< Tangram::Point3 > planet_centers(nplanets);
  double planet_shaft_rad = (sun_shaft_rad/sun_in_rad)*(planet_rad - 0.5*tooth_depth);

  std::vector< Tangram::MatPoly<3> > planet_polys(nplanets), 
    planet_holes(nplanets), planet_shafts(nplanets), planet_shaft_gaskets(nplanets);

  for (int iplanet = 0; iplanet < nplanets; iplanet++) {
    std::vector<Tangram::Point3> planet_pts;

    int itooth = iplanet*sun_nteeth/nplanets;
    Tangram::Vector3 sun2planet = 
      0.5*(sun_pts[4*itooth + 2] + sun_pts[4*itooth + 3]) - gears_center;
    sun2planet.normalize();

    planet_centers[iplanet] = gears_center + 
      (sun_in_rad + gap + planet_rad + 0.5*tooth_depth)*sun2planet;

    planet_pts = cog3d(planet_centers[iplanet], planet_rad - 0.5*tooth_depth, 
      planet_rad + 0.5*tooth_depth, normal, planet_nteeth, dst_tol);

    Tangram::Point3 planet_shaft_center = planet_centers[iplanet] + 
      (carrier_depth + carrier_clearance)*normal;
    std::vector<Tangram::Point3> planet_shaft_pts = 
      circle3d(planet_shaft_center, planet_shaft_rad, normal, 2*shaft_nsides, dst_tol);  

    std::vector<Tangram::Point3> planet_shaft_gasket_pts = 
      circle3d(planet_centers[iplanet], gasket_r_scaling*planet_shaft_rad, 
        normal, 8*shaft_nsides, dst_tol); 

    std::vector<Tangram::Point3> planet_hole_pts = 
      circle3d(planet_centers[iplanet], gasket_r_scaling*planet_shaft_rad, 
        normal, 16*shaft_nsides, dst_tol); 

    Tangram::Vector3 planet_tooth_dir = 
      0.5*(planet_pts[0] + planet_pts[1]) - planet_centers[iplanet];
    planet_tooth_dir.normalize();

    double cos_ang = -Wonton::dot(planet_tooth_dir, sun2planet);
    if (cos_ang > 1.0 - std::numeric_limits<double>::epsilon()) cos_ang = 1.0;
    if (cos_ang < -1.0 + std::numeric_limits<double>::epsilon()) cos_ang = -1.0;
    int sign_rot_angle = Wonton::dot(Wonton::cross(planet_tooth_dir, -sun2planet), 
                                      normal) > 0.0 ? 1 : -1;
    double rot_angle = sign_rot_angle*std::acos(cos_ang);
    Tangram::Matrix rot_matrix = rotation_matrix(normal, rot_angle);

    for (int ivrt = 0; ivrt < planet_pts.size(); ivrt++)
      planet_pts[ivrt] = planet_centers[iplanet] + 
        rot_matrix*(planet_pts[ivrt] - planet_centers[iplanet]);

    planet_polys[iplanet] = prism(planet_pts, gears_depth, 1.0, dst_tol, normal);

    Tangram::Vector3 planet_shaft_dir = planet_shaft_pts[0] - planet_centers[iplanet];
    planet_shaft_dir.normalize();

    cos_ang = -Wonton::dot(planet_shaft_dir, sun2planet);
    if (cos_ang > 1.0 - std::numeric_limits<double>::epsilon()) cos_ang = 1.0;
    if (cos_ang < -1.0 + std::numeric_limits<double>::epsilon()) cos_ang = -1.0;
    sign_rot_angle = Wonton::dot(Wonton::cross(planet_shaft_dir, -sun2planet), 
                                  normal) > 0.0 ? 1 : -1;
    rot_angle = sign_rot_angle*std::acos(cos_ang);
    rot_matrix = rotation_matrix(normal, rot_angle);

    for (int ivrt = 0; ivrt < planet_shaft_pts.size(); ivrt++)
      planet_shaft_pts[ivrt] = planet_shaft_center + 
        rot_matrix*(planet_shaft_pts[ivrt] - planet_shaft_center);  

    planet_shafts[iplanet] = prism(planet_shaft_pts, 
      gears_depth + 2*(carrier_clearance + carrier_depth), 1.0, dst_tol, normal);

    for (int ivrt = 0; ivrt < planet_shaft_gasket_pts.size(); ivrt++)
      planet_shaft_gasket_pts[ivrt] = planet_centers[iplanet] + 
        rot_matrix*(planet_shaft_gasket_pts[ivrt] - planet_centers[iplanet]);  

    planet_shaft_gaskets[iplanet] = prism(planet_shaft_gasket_pts, gears_depth, 1.0, 
                                          dst_tol, normal);

    for (int ivrt = 0; ivrt < planet_hole_pts.size(); ivrt++)
      planet_hole_pts[ivrt] = planet_centers[iplanet] + 
        rot_matrix*(planet_hole_pts[ivrt] - planet_centers[iplanet]);  

    planet_holes[iplanet] = prism(planet_hole_pts, gears_depth, 1.0, dst_tol, normal);
  }

  Tangram::MatPoly<3> carrier_polys[2];
  {
    std::vector< std::vector<Tangram::Point3> > carrier_pts(2);
    for (int icarrier = 0; icarrier < 2; icarrier++)
      carrier_pts[icarrier] = 
        cog3d(carrier_centers[icarrier], sun_out_rad, ring_out_rad, normal, 
              nplanets, dst_tol);

    Tangram::Vector3 new_tooth_dir = planet_centers[0] - gears_center;
    new_tooth_dir.normalize();
    Tangram::Vector3 carrier_tooth_dir = 
      0.5*(carrier_pts[0][0] + carrier_pts[0][1]) - carrier_centers[0];
    carrier_tooth_dir.normalize();

    double cos_ang = Wonton::dot(carrier_tooth_dir, new_tooth_dir);
    if (cos_ang > 1.0 - std::numeric_limits<double>::epsilon()) cos_ang = 1.0;
    if (cos_ang < -1.0 + std::numeric_limits<double>::epsilon()) cos_ang = -1.0;
    int sign_rot_angle = Wonton::dot(Wonton::cross(carrier_tooth_dir, new_tooth_dir), 
                                      normal) > 0.0 ? 1 : -1;
    double rot_angle = sign_rot_angle*std::acos(cos_ang);
    Tangram::Matrix rot_matrix = rotation_matrix(normal, rot_angle);

    for (int icarrier = 0; icarrier < 2; icarrier++)
      for (int ivrt = 0; ivrt < carrier_pts[icarrier].size(); ivrt++)
        carrier_pts[icarrier][ivrt] = carrier_centers[icarrier] + 
          rot_matrix*(carrier_pts[icarrier][ivrt] - carrier_centers[icarrier]); 

    for (int icarrier = 0; icarrier < 2; icarrier++)
      carrier_polys[icarrier] = prism(carrier_pts[icarrier], carrier_depth, 1.0, 
                                      dst_tol, normal);
  }

  Tangram::MatPoly<3> carrier_cen_polys[2];
  {
    std::vector<Tangram::Point3> carrier_cen_pts = 
      circle3d(carrier_centers[0], 0.5*(sun_in_rad + sun_shaft_rad), 
               normal, 4*shaft_nsides, dst_tol);

    Tangram::Vector3 carrier_cen_dir = carrier_cen_pts[0] - carrier_centers[0];
    carrier_cen_dir.normalize();
    Tangram::Vector3 sun2planet = planet_centers[0] - gears_center;
    sun2planet.normalize();

    double cos_ang = -Wonton::dot(carrier_cen_dir, sun2planet);
    if (cos_ang > 1.0 - std::numeric_limits<double>::epsilon()) cos_ang = 1.0;
    if (cos_ang < -1.0 + std::numeric_limits<double>::epsilon()) cos_ang = -1.0;
    double sign_rot_angle = Wonton::dot(Wonton::cross(carrier_cen_dir, -sun2planet), 
                                         normal) > 0.0 ? 1 : -1;
    double rot_angle = sign_rot_angle*std::acos(cos_ang);
    Tangram::Matrix rot_matrix = rotation_matrix(normal, rot_angle);

    for (int ivrt = 0; ivrt < carrier_cen_pts.size(); ivrt++)
      carrier_cen_pts[ivrt] = carrier_centers[0] + 
        rot_matrix*(carrier_cen_pts[ivrt] - carrier_centers[0]);  

    carrier_cen_polys[0] = prism(carrier_cen_pts, carrier_depth, 1.0, dst_tol, normal);

    carrier_cen_pts = circle3d(carrier_centers[1], sun_shaft_rad, normal, 
                               shaft_nsides, dst_tol);
    for (int ivrt = 0; ivrt < carrier_cen_pts.size(); ivrt++)
      carrier_cen_pts[ivrt] = carrier_centers[1] + 
        rot_matrix*(carrier_cen_pts[ivrt] - carrier_centers[1]); 

    carrier_cen_polys[1] = prism(carrier_cen_pts, 
      carrier_depth + carrier_clearance, 1.0, dst_tol, normal);
  }

  Tangram::MatPoly<3> ring_poly;
  {
    std::vector<Tangram::Point3> ring_pts;
    Tangram::Vector3 new_tooth_dir = 0.5*(sun_pts[2] + sun_pts[3]) - gears_center;
    new_tooth_dir.normalize();

    ring_pts = cog3d(gears_center, ring_out_rad - tooth_depth, ring_out_rad, normal, 
                     ring_nteeth, dst_tol);
    Tangram::Vector3 ring_tooth_dir = 0.5*(ring_pts[0] + ring_pts[1]) - gears_center;
    ring_tooth_dir.normalize();

    double cos_ang = Wonton::dot(ring_tooth_dir, new_tooth_dir);
    if (cos_ang > 1.0 - std::numeric_limits<double>::epsilon()) cos_ang = 1.0;
    if (cos_ang < -1.0 + std::numeric_limits<double>::epsilon()) cos_ang = -1.0;
    int sign_rot_angle = Wonton::dot(Wonton::cross(ring_tooth_dir, new_tooth_dir), 
                                      normal) > 0.0 ? 1 : -1;
    double rot_angle = sign_rot_angle*std::acos(cos_ang);
    Tangram::Matrix rot_matrix = rotation_matrix(normal, rot_angle);

    for (int ivrt = 0; ivrt < ring_pts.size(); ivrt++)
      ring_pts[ivrt] = gears_center + rot_matrix*(ring_pts[ivrt] - gears_center);  

    ring_poly = prism(ring_pts, gears_depth, 1.0, dst_tol, normal);
  }

  std::vector<Tangram::Point3> frame_pts = 
    circle3d(frame_center, frame_rad, normal, 4*ring_nteeth, dst_tol);
  Tangram::MatPoly<3> frame_poly = prism(frame_pts, frame_depth, 1.0, dst_tol, normal);

  std::vector<Tangram::Point3> ring_box_pts(4);
  for (int ivrt = 0; ivrt < 4; ivrt++) {
    Tangram::Point3 mid_pt = 0.5*(frame_pts[ivrt*ring_nteeth] + 
      frame_pts[((ivrt + 1)%4)*ring_nteeth]);
    mid_pt += (-2*carrier_clearance - carrier_depth)*normal;
    //Use 1.0e-6 epsilon to make the box slightly wider than the frame
    ring_box_pts[ivrt] = (2.0 + 1.0e-6)*mid_pt + (-1.0 - 1.0e-6)*gears_center;
  }

  Tangram::MatPoly<3> ring_box = prism(ring_box_pts, gears_depth, 1.0, dst_tol, normal);

  res_poly_sets_data.resize(4*nplanets + 9);
  std::vector< std::shared_ptr<RefPolyData_t> > mmpolys_data, rem_polys_data;

  //Cutting out frame interior
  apply_poly(polys_data, frame_poly, mmpolys_data, res_poly_sets_data[0],
             vol_tol, dst_tol, true);

  //Cutting out sun's shaft
  apply_poly(mmpolys_data, sun_shaft_poly, res_poly_sets_data[1], rem_polys_data, 
             vol_tol, dst_tol, true);
  mmpolys_data = rem_polys_data;
  rem_polys_data.clear();

  for (int iplanet = 0; iplanet < nplanets; iplanet++) {
    //Cutting out shaft of planet
    apply_poly(mmpolys_data, planet_shafts[iplanet], 
      res_poly_sets_data[iplanet + 2], rem_polys_data, vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();    
  }    

  for (int iplanet = 0; iplanet < nplanets; iplanet++) {
    //Cutting out shaft gasket of the planet
    apply_poly(mmpolys_data, planet_shaft_gaskets[iplanet], 
      res_poly_sets_data[nplanets + 2 + iplanet], rem_polys_data, vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();    
  } 

  for (int iplanet = 0; iplanet < nplanets; iplanet++) {
    //Cutting out center of planet
    apply_poly(mmpolys_data, planet_holes[iplanet], 
      res_poly_sets_data[2*nplanets + 2 + iplanet], rem_polys_data, vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();    
  } 

  for (int icarrier = 0; icarrier < 2; icarrier++) {
    //Cutting out center of carrier
    apply_poly(mmpolys_data, carrier_cen_polys[icarrier], 
      res_poly_sets_data[3*nplanets + 2 + icarrier], rem_polys_data, vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();    
  }     

  for (int icarrier = 0; icarrier < 2; icarrier++) {
    //Cutting out carrier
    apply_poly(mmpolys_data, carrier_polys[icarrier], 
      res_poly_sets_data[3*nplanets + 4 + icarrier], rem_polys_data, vol_tol, dst_tol, false);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();    
  }      
  std::cout << "Center of the inner carrier is at (" << 
    carrier_centers[0].asV() - 0.5*carrier_depth*normal << ")" << std::endl <<  
    "Center of the outer carrier is at (" << 
    carrier_centers[1].asV() - 0.5*carrier_depth*normal << ")" << std::endl;

  {
    std::vector< std::shared_ptr<RefPolyData_t> > ring_box_poly_data;
    //Cutting out ring box
    apply_poly(mmpolys_data, ring_box, ring_box_poly_data, rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
  
    //Cutting out ring interior from the box
    apply_poly(ring_box_poly_data, ring_poly, rem_polys_data, 
      res_poly_sets_data[3*nplanets + 6], vol_tol, dst_tol, false);

    //Remaining polys are outside the box or inside the ring  
    mmpolys_data.reserve(mmpolys_data.size() + rem_polys_data.size());
    mmpolys_data.insert(mmpolys_data.end(), rem_polys_data.begin(), rem_polys_data.end());
    rem_polys_data.clear();
  }

  //Cutting out the sun gear
  apply_poly(mmpolys_data, sun_poly, res_poly_sets_data[3*nplanets + 7], rem_polys_data, 
             vol_tol, dst_tol, false);
  mmpolys_data = rem_polys_data;
  rem_polys_data.clear();  
  std::cout << "Center of the sun gear is at (" << 
    gears_center.asV() - 0.5*gears_depth*normal << ")" << std::endl;

  for (int iplanet = 0; iplanet < nplanets; iplanet++) {
    //Cutting out planet gear
    apply_poly(mmpolys_data, planet_polys[iplanet], 
      res_poly_sets_data[3*nplanets + 8 + iplanet], rem_polys_data, vol_tol, dst_tol, false);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();    
    std::cout << "Center of the planet gear " << iplanet + 1 << " is at (" << 
      planet_centers[iplanet].asV() - 0.5*gears_depth*normal << ")" << std::endl;
  }
  res_poly_sets_data[4*nplanets + 8] = mmpolys_data;

  return carrier_cen_polys[1];
}

/*!
 @brief For a given mesh computes volume fractions, centroids, and reference material polys
 for the 3D Rotor example. RPGTools and r3d are used to find intersections of reference shapes
 with mesh cells. 
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.

 @param[in] mesh Mesh wrapper.
 @param[out] mesh_material_IDs IDs of materials used in the 3D Rotor example
 @param[out] mesh_material_names Names of materials used in the 3D Rotor example
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
 @param[in] vol_tol Volume tolerance
 @param[in] dst_tol Distance tolerance                                
 @param[in] decompose_cells If mesh has non-convex cells, this flag should be set to true
 in order to decompose cells into tetrahedrons  
 @param[out] reference_mat_polys For every cell and every material inside that cell, 
 a pointer to the collection of single-material polyhedra containing that material
*/
template <class Mesh_Wrapper>
void rotor_material_moments(const Mesh_Wrapper& mesh,
                            std::vector<int>& mesh_material_IDs,
                            std::vector< std::string >& mesh_material_names,
                            std::vector<int>& cell_num_mats,
                            std::vector<int>& cell_mat_ids,
                            std::vector<double>& cell_mat_volfracs,
                            std::vector< Tangram::Point<3> >& cell_mat_centroids,
                            const double vol_tol,
                            const double dst_tol,
                            const bool decompose_cells,
                            std::vector< std::vector< std::vector<r3d_poly> > >*
                              reference_mat_polys = nullptr ) {
  mesh_material_names = {"Water", "Air", "Main shaft", "Blades assembly", "Gears frames",
                         "Ring gears", "Sun gears", "Planet gears", "Carriers", "Outer shafts",
                         "Primary shafts gaskets", "Gear frames sealants", "Planet shafts gaskets"};
  int nmesh_mat_IDs = static_cast<int>(mesh_material_names.size());
  mesh_material_IDs.resize(nmesh_mat_IDs);
  std::iota(mesh_material_IDs.begin(), mesh_material_IDs.end(), 0);                       
   

  std::vector< std::shared_ptr<RefPolyData_t> > mesh_polys;
  mesh_to_r3d_polys<Mesh_Wrapper>(mesh, mesh_polys, dst_tol, decompose_cells);

  std::cout << "Initial number of reference polys: " << mesh_polys.size() << std::endl;

  std::vector<int> sets_material_IDs;
  std::vector< std::vector< std::shared_ptr<RefPolyData_t> > > ref_poly_sets;

  double chamber_in_rad = 0.25;
  double chamber_out_rad = 0.26;
  double dchamber_r = chamber_out_rad - chamber_in_rad;
  Tangram::Point3 shaft_centroid(0.5, 0.5, 0.5);
  double shaft_rad = 0.02;
  int nshaft_sides = 6;
  double shaft_len = 2.5*chamber_out_rad;
  double carrier_shaft_len = 0.25*chamber_out_rad;
  double gasket_r_scaling = 1.02;
  
  double shaft_mount_rad = 1.25*shaft_rad;
  int nshaft_mount_sides = 12;
  double blades_mount_rad = 1.5*shaft_rad;
  int nblades = 5;
  double blades_mount_len = 0.6*chamber_in_rad;
  double blades_height = 0.95*chamber_in_rad - blades_mount_rad;
  double blades_base_scaling = blades_mount_rad/(blades_height + blades_mount_rad);

  double sun_shaft_rad[2] = {0.75*shaft_rad, 0.75*shaft_rad};
  double sun_in_rad[2] = {shaft_rad, 0.95*shaft_rad};
  double sun_out_rad[2] = {1.25*shaft_rad, 1.06875*shaft_rad};
  int sun_nteeth[2] = {15, 18};
  int nplanets[2] = {5, 3};
  int planet_nteeth[2] = {6, 10};
  double gear_ring_rad[2] = {2*sun_out_rad[0], 2.4*sun_out_rad[1]};
  double gear_frame_rad[2] = {gear_ring_rad[0] + dchamber_r, gear_ring_rad[1] + dchamber_r};
  int ngear_frame_sides = 12;

  Tangram::Point3 gear_frame_joint_base_cen[2] = {shaft_centroid, shaft_centroid};
  gear_frame_joint_base_cen[0][1] -= sqrt(pow2(chamber_in_rad) - 0.25*pow2(gear_frame_rad[0]));
  gear_frame_joint_base_cen[1][1] += sqrt(pow2(chamber_in_rad) - 0.25*pow2(gear_frame_rad[1]));
  double gear_frame_joint_height[2] = {
    chamber_out_rad - (shaft_centroid[1] - gear_frame_joint_base_cen[0][1]), 
    chamber_out_rad - (gear_frame_joint_base_cen[1][1] - shaft_centroid[1])};

  double gear_frame_outer_height[2] = {0.5*dchamber_r, 0.5*dchamber_r};
  Tangram::Point3 gear_frame_outer_base_cen[2] = {shaft_centroid, shaft_centroid};
  gear_frame_outer_base_cen[0][1] -= 0.5*shaft_len - gear_frame_outer_height[0];
  gear_frame_outer_base_cen[1][1] += 0.5*shaft_len - gear_frame_outer_height[1];

  Tangram::Point3 gear_frame_main_base_cen[2] = {shaft_centroid, shaft_centroid};
  gear_frame_main_base_cen[0][1] -= chamber_out_rad;
  gear_frame_main_base_cen[1][1] += chamber_out_rad;
  double gear_frame_main_height[2] = {
    0.5*shaft_len - gear_frame_outer_height[0] - (shaft_centroid[1] - gear_frame_main_base_cen[0][1]),
    0.5*shaft_len - gear_frame_outer_height[1] - (gear_frame_main_base_cen[1][1] - shaft_centroid[1])};
  double gear_carrier_clearance[2] = {
    0.075*gear_frame_main_height[0], 0.075*gear_frame_main_height[1]};
  double gear_carrier_depth[2] = {
    0.25*(gear_frame_main_height[0] - 4*gear_carrier_clearance[0]),
    0.25*(gear_frame_main_height[1] - 4*gear_carrier_clearance[1])};

  Tangram::Vector3 shaft_base_normal(0.0, 1.0, 0.0);

  //Temporary sets
  std::vector< std::shared_ptr<RefPolyData_t> > mmpolys_data, rem_polys_data;

  //Prepare shaft's central part
  double shaft_cen_len = 2*(shaft_centroid[1] - gear_frame_main_base_cen[0][1]);
  std::vector<Tangram::Point3> shaft_cen_base = 
    circle3d(shaft_centroid + Tangram::Point3(0.0, 0.5*shaft_cen_len, 0.0), 
             shaft_rad, shaft_base_normal, nshaft_sides, dst_tol);
  Tangram::MatPoly<3> shaft_cen_poly = 
    prism(shaft_cen_base, shaft_cen_len, 1.0, dst_tol, shaft_base_normal);             

  std::cout << "Preparing the central part..." << std::endl;
  int icur_ref_set = ref_poly_sets.size();
  ref_poly_sets.resize(icur_ref_set + 1);

  apply_poly(mesh_polys, shaft_cen_poly, ref_poly_sets[icur_ref_set], mmpolys_data, 
             vol_tol, dst_tol, true);
  sets_material_IDs.push_back(mesh_material_IDs[2]);

  mesh_polys.clear();

  //Prepare shaft's mounting part
  std::vector<Tangram::Point3> shaft_mount_base = 
    circle3d(shaft_centroid + Tangram::Point3(0.0, 0.5*blades_mount_len, 0.0), 
             shaft_mount_rad, shaft_base_normal, nshaft_mount_sides, dst_tol);
  Tangram::MatPoly<3> shaft_mount_poly = 
    prism(shaft_mount_base, blades_mount_len, 1.0, dst_tol, shaft_base_normal);   

  icur_ref_set = ref_poly_sets.size();
  ref_poly_sets.resize(icur_ref_set + 1);

  apply_poly(mmpolys_data, shaft_mount_poly, ref_poly_sets[icur_ref_set], rem_polys_data, 
             vol_tol, dst_tol, true);
  mmpolys_data = rem_polys_data;
  rem_polys_data.clear();
  sets_material_IDs.push_back(mesh_material_IDs[2]);

  //Prepare blades' mounting part
  std::vector<Tangram::Point3> blades_mount_base = 
    circle3d(shaft_centroid + Tangram::Point3(0.0, 0.5*blades_mount_len, 0.0), 
             blades_mount_rad, shaft_base_normal, 2*nblades, dst_tol);
  Tangram::MatPoly<3> blades_mount_poly = 
    prism(blades_mount_base, blades_mount_len, 1.0, dst_tol, shaft_base_normal);   

  icur_ref_set = ref_poly_sets.size();
  ref_poly_sets.resize(icur_ref_set + 1);

  apply_poly(mmpolys_data, blades_mount_poly, ref_poly_sets[icur_ref_set], rem_polys_data,
             vol_tol, dst_tol, true);
  mmpolys_data = rem_polys_data;
  rem_polys_data.clear();
  sets_material_IDs.push_back(mesh_material_IDs[3]);

  for (int iblade = 0; iblade < nblades; iblade++) {
    //Prepare blade
    std::vector<Tangram::Point3> blade_base;
    blade_base.reserve(4);
    const std::vector<int>& mount_vrts = blades_mount_poly.face_vertices(2*iblade + 1);
    assert(mount_vrts.size() == 4);
    for (int ivrt = 0; ivrt < 4; ivrt++)
      blade_base.push_back(blades_mount_poly.vertex_point(mount_vrts[3 - ivrt]));

    Tangram::MatPoly<3> blade_poly = 
      prism(blade_base, blades_height, blades_base_scaling, dst_tol);   

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, blade_poly, ref_poly_sets[icur_ref_set], rem_polys_data, 
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[3]);
  }  

  std::cout << "Center of blades assembly is at (" << shaft_centroid << ")" << std::endl;

  for (int igear = 0; igear < 2; igear++) {
    std::cout << "Preparing gearset " << igear + 1 << "..." << std::endl;
    Tangram::Vector3 gears_base_normal = (igear == 0) ? shaft_base_normal : -shaft_base_normal;

    //Prepare main shaft gasket
    std::vector<Tangram::Point3> shaft_gasket_base = 
      circle3d(gear_frame_joint_base_cen[igear], gasket_r_scaling*shaft_rad, 
        gears_base_normal, 4*nshaft_sides, dst_tol);
    Tangram::MatPoly<3> shaft_gasket_poly = 
      prism(shaft_gasket_base, gear_frame_joint_height[igear], 1.0, 
            dst_tol, gears_base_normal);

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, shaft_gasket_poly, ref_poly_sets[icur_ref_set], rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[10]);

    //Prepare joint part, sealant
    std::vector<Tangram::Point3> gear_frame_joint_base = 
      circle3d(gear_frame_joint_base_cen[igear], (2*gasket_r_scaling - 1.0)*shaft_rad, 
        gears_base_normal, 8*nshaft_sides, dst_tol);
    Tangram::MatPoly<3> gear_frame_joint_poly = 
      prism(gear_frame_joint_base, gear_frame_joint_height[igear], 1.0, 
            dst_tol, gears_base_normal);   

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, gear_frame_joint_poly, ref_poly_sets[icur_ref_set], rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[11]);

    //Prepare joint part, outer frame
    gear_frame_joint_base = circle3d(gear_frame_joint_base_cen[igear], 
      gear_frame_rad[igear], gears_base_normal, ngear_frame_sides, dst_tol);
    gear_frame_joint_poly = prism(gear_frame_joint_base, 
      gear_frame_joint_height[igear], 1.0, dst_tol, gears_base_normal);

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, gear_frame_joint_poly, ref_poly_sets[icur_ref_set], rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[4]);

    std::cout << "Center of the inner part of the frame is at (" << 
      gear_frame_joint_base_cen[igear].asV() - 
      0.5*gear_frame_joint_height[igear]*gears_base_normal << ")" << std::endl;

    //Prepare gears part, outer frame
    std::vector<Tangram::Point3> gear_frame_main_base = circle3d(gear_frame_main_base_cen[igear], 
      gear_frame_rad[igear], gears_base_normal, ngear_frame_sides, dst_tol);
    Tangram::MatPoly<3> gear_frame_main_poly = prism(gear_frame_main_base, 
      gear_frame_main_height[igear], 1.0, dst_tol, gears_base_normal); 

    std::vector< std::shared_ptr<RefPolyData_t> > planetary_gear_polys;

    apply_poly(mmpolys_data, gear_frame_main_poly, planetary_gear_polys, rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();

    std::vector< std::vector< std::shared_ptr<RefPolyData_t> > > planetary_gear_poly_sets;
    Tangram::MatPoly<3> gear_carrier_shaft_poly = 
      planetary_gear(planetary_gear_polys, gears_base_normal, gear_frame_main_base_cen[igear],
        gear_frame_main_height[igear], sun_shaft_rad[igear], sun_in_rad[igear], sun_out_rad[igear], 
        gear_carrier_clearance[igear], gear_carrier_depth[igear], gear_ring_rad[igear], 
        nplanets[igear], sun_nteeth[igear], planet_nteeth[igear], nshaft_sides, gasket_r_scaling,
        vol_tol, dst_tol, planetary_gear_poly_sets);

    //Gearset outer frame
    ref_poly_sets.push_back(planetary_gear_poly_sets[0]);
    sets_material_IDs.push_back(mesh_material_IDs[4]);

    //Gearset sun's shaft
    ref_poly_sets.push_back(planetary_gear_poly_sets[1]);
    sets_material_IDs.push_back(mesh_material_IDs[2]);

    //Gearset planets' shafts
    for (int iplanet = 0; iplanet < nplanets[igear]; iplanet++) {
      ref_poly_sets.push_back(planetary_gear_poly_sets[iplanet + 2]);
      sets_material_IDs.push_back(mesh_material_IDs[8]);
    }

    //Gearset planets' shafts gaskets
    for (int iplanet = 0; iplanet < nplanets[igear]; iplanet++) {
      ref_poly_sets.push_back(planetary_gear_poly_sets[nplanets[igear] + 2 + iplanet]);
      sets_material_IDs.push_back(mesh_material_IDs[12]);
    }

    //Gearset planets' centers
    for (int iplanet = 0; iplanet < nplanets[igear]; iplanet++) {
      ref_poly_sets.push_back(planetary_gear_poly_sets[2*nplanets[igear] + 2 + iplanet]);
      sets_material_IDs.push_back(mesh_material_IDs[1]);
    }    

    //Gearset inner carrier center
    ref_poly_sets.push_back(planetary_gear_poly_sets[3*nplanets[igear] + 2]);
    sets_material_IDs.push_back(mesh_material_IDs[1]);

    //Gearset outer carrier shaft
    ref_poly_sets.push_back(planetary_gear_poly_sets[3*nplanets[igear] + 3]);
    sets_material_IDs.push_back(mesh_material_IDs[9]);

    //Gearset carriers
    for (int icarrier = 0; icarrier < 2; icarrier++) {
      ref_poly_sets.push_back(planetary_gear_poly_sets[3*nplanets[igear] + 4 + icarrier]);
      sets_material_IDs.push_back(mesh_material_IDs[8]);      
    }

    //Gearset ring
    ref_poly_sets.push_back(planetary_gear_poly_sets[3*nplanets[igear] + 6]);
    sets_material_IDs.push_back(mesh_material_IDs[5]);

    //Gearset sun
    ref_poly_sets.push_back(planetary_gear_poly_sets[3*nplanets[igear] + 7]);
    sets_material_IDs.push_back(mesh_material_IDs[6]);

    //Gearset planets
    for (int iplanet = 0; iplanet < nplanets[igear]; iplanet++) {
      ref_poly_sets.push_back(planetary_gear_poly_sets[3*nplanets[igear] + 8 + iplanet]);
      sets_material_IDs.push_back(mesh_material_IDs[7]);
    }

    //Gearset gaps
    ref_poly_sets.push_back(planetary_gear_poly_sets[4*nplanets[igear] + 8]);
    sets_material_IDs.push_back(mesh_material_IDs[1]);

    //Prepare carrier shaft
    std::vector< Tangram::Point<3> > carrier_shaft_pts = 
      gear_carrier_shaft_poly.face_points(nshaft_sides + 1);
    assert(carrier_shaft_pts.size() == nshaft_sides);
    std::reverse(carrier_shaft_pts.begin(), carrier_shaft_pts.end());
 
    Tangram::MatPoly<3> carrier_shaft_poly = 
      prism(carrier_shaft_pts, carrier_shaft_len, 1.0, dst_tol, gears_base_normal);

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, carrier_shaft_poly, ref_poly_sets[icur_ref_set], rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[9]);

    //Prepare carrier shaft gasket
    std::vector<Tangram::Point3> carrier_shaft_gasket_base = 
      circle3d(gear_frame_outer_base_cen[igear], gasket_r_scaling*sun_shaft_rad[igear], 
        gears_base_normal, 4*nshaft_sides, dst_tol);
    Tangram::MatPoly<3> carrier_shaft_gasket_poly = 
      prism(carrier_shaft_gasket_base, gear_frame_outer_height[igear], 1.0, 
            dst_tol, gears_base_normal);

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, carrier_shaft_gasket_poly, 
      ref_poly_sets[icur_ref_set], rem_polys_data, vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[10]);

    //Prepare outer frame sealant
    std::vector<Tangram::Point3> gear_frame_outer_base = 
      circle3d(gear_frame_outer_base_cen[igear], (2*gasket_r_scaling - 1.0)*sun_shaft_rad[igear], 
        gears_base_normal, 8*nshaft_sides, dst_tol);
    Tangram::MatPoly<3> gear_frame_outer_poly = 
      prism(gear_frame_outer_base, gear_frame_outer_height[igear], 1.0, 
            dst_tol, gears_base_normal);   

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, gear_frame_outer_poly, ref_poly_sets[icur_ref_set], rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[11]);

    //Prepare outer frame
    gear_frame_outer_base = circle3d(gear_frame_outer_base_cen[igear], 
      gear_frame_rad[igear], gears_base_normal, ngear_frame_sides, dst_tol);
    gear_frame_outer_poly = prism(gear_frame_outer_base, 
      gear_frame_outer_height[igear], 1.0, dst_tol, gears_base_normal);   

    icur_ref_set = ref_poly_sets.size();
    ref_poly_sets.resize(icur_ref_set + 1);

    apply_poly(mmpolys_data, gear_frame_outer_poly, ref_poly_sets[icur_ref_set], rem_polys_data,
               vol_tol, dst_tol, true);
    mmpolys_data = rem_polys_data;
    rem_polys_data.clear();
    sets_material_IDs.push_back(mesh_material_IDs[4]);

    std::cout << "Center of the outer part of the frame is at (" << 
      gear_frame_outer_base_cen[igear].asV() - 
      0.5*gear_frame_outer_height[igear]*gears_base_normal << ")" << std::endl;
  }

  //Making all exterior polys set
  ref_poly_sets.push_back(mmpolys_data);
  sets_material_IDs.push_back(mesh_material_IDs[0]);

  std::cout << "Finalizing rotor data..." << std::endl;
  finalize_ref_data(mesh, ref_poly_sets, sets_material_IDs, cell_num_mats, cell_mat_ids,
    cell_mat_volfracs, cell_mat_centroids, dst_tol, decompose_cells, reference_mat_polys);
  std::cout << "Done with rotor data!" << std::endl;
}

#endif
