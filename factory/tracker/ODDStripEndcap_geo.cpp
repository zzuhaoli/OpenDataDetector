
// Open Data Dector project
//
// (c) 2021 CERN for the benefit of the ODD project
//
// Mozilla Public License Version 2.0

#ifdef ODD_ACTS_EXTENSION
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#endif

#include "DD4hep/DetFactoryHelper.h"
#include "ODDModuleHelper.hpp"
#include "ODDServiceHelper.hpp"

#include <vector>

using namespace std;
using namespace dd4hep;

/// Standard create_element(...) DD4hep method for Strip endcaps
///
/// This can build long/short strip endcaps for the Open Data Detector,
/// allowing a variable number of endcap discs that consist of a
/// variable number of rings each.
///
/// It supports double (with stereo angle) and single sided modules,
/// the latter are assumed to be trapezoidal.
///
/// @param oddd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens the sensitive detector descrition
///
/// @return a reference counted DetElement
static Ref_t create_element(Detector &oddd, xml_h xml, SensitiveDetector sens)
{
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make DetElement
  DetElement endcapDetector(detName, x_det.id());

#ifdef ODD_ACTS_EXTENSION
  // Add Extension to DetElement for the RecoGeometry
  Acts::ActsExtension *endcapExtension = new Acts::ActsExtension();
  endcapExtension->addType("endcap", "detector");

  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat)
  {
    xml_comp_t x_boundary_material = bmat;
    xmlToProtoSurfaceMaterial(x_boundary_material, *endcapExtension,
                              "boundary_material");
  }
  endcapDetector.addExtension<Acts::ActsExtension>(endcapExtension);
#endif

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  string endcapShapeName = x_det_dim.nameStr();

  Tube endcapShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume endcapVolume(detName, endcapShape, oddd.air());
  endcapVolume.setVisAttributes(oddd, x_det.visStr());

  Assembly diskAssembly("Disk");

  // DetElement tree
  DetElement diskElementTemplate("DiskElementTemplate", 0);

  // build the ring templates
  size_t ringNum = 0;
  for (xml_coll_t ring(x_det, _U(ring)); ring; ++ring, ++ringNum)
  {
    xml_comp_t x_ring = ring;

    string ringName = "Ring" + std::to_string(ringNum);
    Assembly ringAssembly(ringName);

    // DetElement tree
    DetElement ringElement(ringName, ringNum);

    // Build the module
    if (x_ring.hasChild(_U(module)))
    {
      xml_comp_t x_module = x_ring.child(_U(module));
      auto module =
          ODDModuleHelper::assembleTrapezoidalModule(oddd, sens, x_module);

      // Place the modules
      unsigned int nPhi = x_ring.nphi();
      double phiStep = 2 * M_PI / nPhi;
      double phi0 = x_ring.phi0();
      double r = x_ring.r();
      double zgap = x_ring.gap();

      for (unsigned int modNum = 0; modNum < nPhi; ++modNum)
      {
        // The module name
        string moduleName = _toString((int)modNum, "module%d");

        bool odd = bool(modNum % 2);

        // Position parameters
        double phi = phi0 + modNum * phiStep;
        double x = r * cos(phi);
        double y = r * sin(phi);
        double z = odd ? -zgap : zgap;

        // Place Module Box Volumes, flip if necessary
        Position trans(x, y, z);
        double flip = odd ? M_PI : 0.;

        double angX = 0.5 * M_PI + flip;
        double angY = odd ? 0.5 * M_PI - phi : 0.5 * M_PI + phi;

        PlacedVolume placedModule = ringAssembly.placeVolume(
            module.first,
            Transform3D(
                RotationX(angX) * RotationY(angY),
                trans)); // RotationZ(phi + 1.5 * M_PI) * RotationY(flip)
        placedModule.addPhysVolID("module", modNum);
        // Clone the detector element
        auto moduleElement = module.second.clone(moduleName, modNum);
        moduleElement.setPlacement(placedModule);
        // Assign it as child to the stave template
        ringElement.add(moduleElement);
      }

      // Now add the ring detector Element to the disk
      diskElementTemplate.add(ringElement);

      size_t supportNum = 0;
      for (xml_coll_t sup(x_ring, _U(support)); sup; ++sup, ++supportNum)
      {
        xml_comp_t x_support = sup;
        // Create the volume of the support structure
        string supportName = _toString((int)supportNum, "RingSupport%d");
        Volume supportVolume(
            supportName,
            Tube(x_support.rmin(), x_support.rmax(), x_support.dz()),
            oddd.material(x_support.materialStr()));
        supportVolume.setVisAttributes(oddd, x_support.visStr());
        // Place the support structure
        ringAssembly.placeVolume(supportVolume,
                                 Position(0., 0., x_support.z_offset()));
      }

      // Cooling rings
      buildCoolingRings(oddd, ringAssembly, x_ring);

      PlacedVolume placedRing = diskAssembly.placeVolume(
          ringAssembly, Position(0., 0., x_ring.z_offset()));

      placedRing.addPhysVolID("ring", ringNum);
      ringElement.setPlacement(placedRing);
    }
  }

  // Loop over the layers and place the disk, remember the z positions
  std::vector<double> endcapZ;
  size_t layNum = 0;
  for (xml_coll_t lay(xml, _U(layer)); lay; ++lay, ++layNum)
  {
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
    PlacedVolume placedLayer =
        endcapVolume.placeVolume(layerVolume, Position(0., 0., zeff));
    placedLayer.addPhysVolID("layer", layNum);
    endcapZ.push_back(zeff);

#ifdef ODD_ACTS_EXTENSION
    // Place the layer with appropriate Acts::Extension
    // Configure the ACTS extension
    Acts::ActsExtension *layerExtension = new Acts::ActsExtension();
    layerExtension->addType("sensitive disk", "layer");
    layerExtension->addValue(25, "r_min", "envelope");
    layerExtension->addValue(10, "r_max", "envelope");
    layerExtension->addValue(10, "z_min", "envelope");
    layerExtension->addValue(10, "z_max", "envelope");

    // Check if the disk has a surface binning instruction
    if (x_layer.hasChild(_Unicode(surface_binning)))
    {
      xml_comp_t sfBinning = x_layer.child(_Unicode(surface_binning));
      layerExtension->addValue(sfBinning.attr<int>("nr"), "n_r",
                               "surface_binning");
      layerExtension->addValue(sfBinning.attr<int>("nphi"), "n_phi",
                               "surface_binning");
    }

    // Add the proto layer material
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat)
    {
      xml_comp_t x_layer_material = lmat;
      xmlToProtoSurfaceMaterial(x_layer_material, *layerExtension,
                                "layer_material");
    }
    layerElement.addExtension<Acts::ActsExtension>(layerExtension);

#endif

    // Finish up the DetElement tree
    layerElement.setPlacement(placedLayer);
    endcapDetector.add(layerElement);
  }

  if (x_det.hasChild(_Unicode(services)))
  {
    // Grab the services - cables
    xml_comp_t x_services = x_det.child(_Unicode(services));
    for (xml_coll_t crout(x_services, _Unicode(cable_routing)); crout;
         ++crout)
    {
      xml_comp_t x_cable_routing = crout;
      buildEndcapRouting(oddd, endcapVolume, x_cable_routing, endcapZ);
    }
    // Grab for services - cooling
    for (xml_coll_t crout(x_services, _Unicode(cooling_routing)); crout;
         ++crout)
    {
      xml_comp_t x_cooling_routing = crout;
      buildEndcapRouting(oddd, endcapVolume, x_cooling_routing, endcapZ);
    }
  }

  // Place Volume
  Volume motherVolume = oddd.pickMotherVolume(endcapDetector);
  Position translation(0., 0., x_det_dim.z());
  PlacedVolume placedEndcap =
      motherVolume.placeVolume(endcapVolume, translation);
  // "system" is hard coded in the DD4Hep::VolumeManager
  placedEndcap.addPhysVolID("system", endcapDetector.id());
  endcapDetector.setPlacement(placedEndcap);

  // And return it
  return endcapDetector;
}

DECLARE_DETELEMENT(ODDStripEndcap, create_element)
