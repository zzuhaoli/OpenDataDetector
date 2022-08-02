// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#include <vector>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "ODDHelper.hpp"
#include "ODDModuleHelper.hpp"
#include "ODDServiceHelper.hpp"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;

/// Standard create_element(...) DD4hep method for Pixel endcap
///
/// @param oddd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens the sensitive detector descrition
///
/// @return a reference counted DetElement
static Ref_t create_element(Detector &oddd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make DetElement
  DetElement endcapDetector(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, endcapDetector);

  auto &params = ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
      endcapDetector);

  // Add Extension to DetElement for the RecoGeometry
  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    ODDHelper::xmlToProtoSurfaceMaterial(x_boundary_material, params,
                                         "boundary_material");
  }

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  string endcapShapeName = x_det_dim.nameStr();

  Tube endcapShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume endcapVolume(detName, endcapShape, oddd.air());
  endcapVolume.setVisAttributes(oddd, x_det.visStr());

  // Place Volume
  Volume motherVolume = oddd.pickMotherVolume(endcapDetector);
  Position translation(0., 0., x_det_dim.z());
  PlacedVolume placedEndcap =
      motherVolume.placeVolume(endcapVolume, translation);

  Assembly diskAssembly("disk");

  // DetElement tree
  DetElement diskElementTemplate("DiskElementTemplate", 0);

  // Loop over the rings to create a template disk
  size_t ringNum = 0;
  for (xml_coll_t ring(xml, _U(ring)); ring; ++ring, ++ringNum) {
    // Get the ring
    xml_comp_t x_ring = ring;

    // The ring name
    string ringName = _toString((int)ringNum, "ring%d");
    Assembly ringAssembly(ringName);
    ringAssembly.setVisAttributes(oddd, x_ring.visStr());

    // DetElement tree
    DetElement ringElement(ringName, ringNum);

    if (x_ring.hasChild(_U(module))) {
      xml_comp_t x_module = x_ring.child(_U(module));
      auto module =
          ODDModuleHelper::assembleTrapezoidalModule(oddd, sens, x_module);

      double r = x_ring.r();
      double phi0 = x_ring.phi0();
      unsigned int nModules = x_ring.nphi();
      double zgap = x_ring.gap();
      double phiStep = 2. * M_PI / nModules;

      // Loop over modules
      for (unsigned int modNum = 0; modNum < nModules; ++modNum) {
        // The module name
        string moduleName = _toString((int)modNum, "module%d");

        bool odd = bool(modNum % 2);

        double phi = phi0 + modNum * phiStep;
        double x = r * cos(phi);
        double y = r * sin(phi);
        double z = odd ? -zgap : zgap;

        // Place Module Box Volumes, flip if necessary
        Position trans(x, y, z);

        double angX = 1.5 * M_PI;
        double angY = 0.5 * M_PI - phi;

        PlacedVolume placedModule = ringAssembly.placeVolume(
            module.first,
            Transform3D(RotationX(angX) * RotationY(angY), trans));
        placedModule.addPhysVolID("module", modNum);
        // Clone the detector element
        auto moduleElement = module.second.clone(moduleName, modNum);
        moduleElement.setPlacement(placedModule);
        // Assign it as child to the stave template
        ringElement.add(moduleElement);
      }
    }

    // Place Ring assembly into disk
    PlacedVolume placedRing = diskAssembly.placeVolume(
        ringAssembly, Position(0., 0., x_ring.z_offset()));
    placedRing.addPhysVolID("ring", ringNum);
    ringElement.setPlacement(placedRing);
    // Add it to the Disk element template
    diskElementTemplate.add(ringElement);
  }

  xml_comp_t x_support = x_det.child(_Unicode(ring_support));
  // The support shape
  Tube supportShape(x_support.rmin(), x_support.rmax(), x_support.dz());
  Volume supportVolume("DiskSupport", supportShape,
                       oddd.material(x_support.materialStr()));
  supportVolume.setVisAttributes(oddd, x_support.visStr());
  diskAssembly.placeVolume(supportVolume);

  // Cooling rings
  buildCoolingRings(oddd, diskAssembly, x_det);

  // Loop over the layers and place the disk
  size_t layNum = 0;
  // Remember the layers for the service routing
  std::vector<double> endcapZ;
  for (xml_coll_t lay(xml, _U(layer)); lay; ++lay, ++layNum) {
    // Get the layer
    xml_comp_t x_layer = lay;

    // The Layer envelope volume
    string layerName = detName + std::to_string(layNum);
    Volume layerVolume(layerName,
                       Tube(x_layer.rmin(), x_layer.rmax(), x_layer.dz()),
                       oddd.air());

    layerVolume.setVisAttributes(oddd, x_layer.visStr());

    string diskElName = _toString((int)layNum, "disk%d");

    // The DetElement tree
    DetElement layerElement(layerName, layNum);
    auto diskElement = diskElementTemplate.clone(diskElName, layNum);

    // Place the disk into the layer
    PlacedVolume placedDisk = layerVolume.placeVolume(diskAssembly);
    diskElement.setPlacement(placedDisk);
    layerElement.add(diskElement);

    // Place Ring assembly into disk
    double zeff = x_layer.z_offset() - x_det_dim.z();
    endcapZ.push_back(zeff);

    PlacedVolume placedLayer =
        endcapVolume.placeVolume(layerVolume, Position(0., 0., zeff));
    placedLayer.addPhysVolID("layer", layNum);

    auto &layerParams =
        ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
            layerElement);

    layerParams.set<double>("envelope_r_min", 10.0);
    layerParams.set<double>("envelope_r_max", 10.0);
    layerParams.set<double>("envelope_z_min", 10.0);
    layerParams.set<double>("envelope_z_max", 10.0);

    // Check if the disk has a surface binning instruction
    if (x_layer.hasChild(_Unicode(surface_binning))) {
      xml_comp_t sfBinning = x_layer.child(_Unicode(surface_binning));
      layerParams.set<bool>("surface_binning", true);
      layerParams.set<int>("surface_binning_n_r", sfBinning.attr<int>("nr"));
      layerParams.set<int>("surface_binning_n_phi",
                           sfBinning.attr<int>("nphi"));
    }

    // Add the proto layer material
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      ODDHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                           "layer_material");
    }

    // Finish up the DetElement tree
    layerElement.setPlacement(placedLayer);
    endcapDetector.add(layerElement);
  }

  // Close up the detector
  if (x_det.hasChild(_U(disk))) {
    // Endplate disk
    xml_comp_t x_endplate = x_det.child(_U(disk));

    // The Shape and Volume
    Tube endplateShape(x_endplate.rmin(), x_endplate.rmax(), x_endplate.dz());
    Volume endplateVolume("Endplate", endplateShape,
                          oddd.material(x_endplate.materialStr()));
    endplateVolume.setVisAttributes(oddd, x_endplate.visStr());

    double zeff = x_endplate.z_offset() - x_det_dim.z();
    endcapZ.push_back(zeff);
    PlacedVolume placedEndplate =
        endcapVolume.placeVolume(endplateVolume, Position(0., 0., zeff));

    DetElement endplateElement(x_endplate.nameStr(), 0);
    dd4hep::DetType typeFlags{};
    endplateElement.setTypeFlag(typeFlags.to_ulong());

    // Place the layer with appropriate Acts::Extension
    // Configure the ACTS extension
    auto &layerParams =
        ODDHelper::ensureExtension<dd4hep::rec::VariantParameters>(
            endplateElement);
    // Add the proto layer material
    for (xml_coll_t lmat(x_endplate, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      ODDHelper::xmlToProtoSurfaceMaterial(x_layer_material, layerParams,
                                           "layer_material");
    }

    // Finish up the DetElement tree
    endplateElement.setPlacement(placedEndplate);
    endcapDetector.add(endplateElement);
  }

  if (x_det.hasChild(_Unicode(services))) {
    // Grab the services
    xml_comp_t x_services = x_det.child(_Unicode(services));
    if (x_services.hasChild(_Unicode(cable_routing))) {
      xml_comp_t x_cable_routing = x_services.child(_Unicode(cable_routing));
      buildEndcapRouting(oddd, endcapVolume, x_cable_routing, endcapZ);
    }
    if (x_services.hasChild(_Unicode(cooling_routing))) {
      xml_comp_t x_cooling_routing =
          x_services.child(_Unicode(cooling_routing));
      buildEndcapRouting(oddd, endcapVolume, x_cooling_routing, endcapZ);
    }
  }

  // Place the additional support cylinders per detector
  std::vector<double> layerR;
  buildSupportCylinder(oddd, endcapVolume, x_det, layerR);

  // "system" is hard coded in the DD4Hep::VolumeManager
  placedEndcap.addPhysVolID("system", endcapDetector.id());
  endcapDetector.setPlacement(placedEndcap);

  // And return it
  return endcapDetector;
}

DECLARE_DETELEMENT(ODDPixelEndcap, create_element)
