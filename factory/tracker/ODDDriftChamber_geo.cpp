#include <vector>
#include <iostream>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"

using namespace dd4hep;

struct DCSuperLaylerCalculator {
    // INPUT
    int superlayer_id;
    std::string type; // A/U/V
    int n; // number of layers in the superlayer
    double rmin; // rmin of the superlayer
    double rmax; // rmax of the superlayer

    // OUTPUT

    DCSuperLaylerCalculator(int _id, const std::string& _type, int _n, double _rmin, double _rmax)
    : superlayer_id(_id), type(_type), n(_n), rmin(_rmin), rmax(_rmax) {

    }

    void calc() {
        // According to the rmin/rmax and n, calculate the size of cell (inner)
        double cell_size = (rmax-rmin) / n;

    }

    void dump() {
        std::cout 
        << "SUPER LAYER [" << superlayer_id << "]"
        << " type: " << type
        << " n: " << n
        << " rmin: " << rmin
        << " rmax: " << rmax
        << std::endl;
    }

};

static Ref_t create_element(Detector &oddd, xml_h xml,
                                SensitiveDetector sens) {

    xml_det_t x_det = xml;
    std::string detName = x_det.nameStr();

    // Make Volume
    xml_comp_t x_det_tubs = x_det.child(_U(tubs));

    // Make DetElement
    DetElement driftChamberElement(detName, x_det.id());

    // Envelope
    std::string shapeName = x_det_tubs.nameStr();
    Tube tubeShape(shapeName, x_det_tubs.rmin(), x_det_tubs.rmax(),
                    x_det_tubs.dz());
    Volume tubeVolume(detName, tubeShape,
                        oddd.material(x_det_tubs.materialStr()));
    tubeVolume.setVisAttributes(oddd, x_det.visStr());

    // Place it in the mother
    Volume motherVolume = oddd.pickMotherVolume(driftChamberElement);
    PlacedVolume placedTube = motherVolume.placeVolume(tubeVolume);
    placedTube.addPhysVolID(detName, driftChamberElement.id());
    driftChamberElement.setPlacement(placedTube);

    // Construct super layers
    std::cout << "START BUILDING SUPER LAYER" << std::endl;
    size_t superlayerNum = 0;
    if (x_det.hasChild(_U(layers))) {
        for (xml_coll_t superlay(xml, _U(layers)); 
            superlay; ++superlay, ++superlayerNum) {
            xml_comp_t x_superlayer = superlay;

            DCSuperLaylerCalculator dc_superlayer_calc(superlayerNum,
                                                       x_superlayer.typeStr(),
                                                       x_superlayer.number(),
                                                       x_superlayer.rmin(),
                                                       x_superlayer.rmax());
            dc_superlayer_calc.dump();
        }
    }


  return driftChamberElement;

}

DECLARE_DETELEMENT(ODDDriftChamber, create_element)
