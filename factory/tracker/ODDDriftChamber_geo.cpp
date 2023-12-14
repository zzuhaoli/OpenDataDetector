#include <vector>
#include <tuple>
#include <iostream>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"

using namespace dd4hep;

struct DCSuperLaylerCalculator {
    // INPUT
    const int superlayer_id;
    const std::string type; // A/U/V
    const int n; // number of layers in the superlayer
    const double rmin; // rmin of the superlayer
    const double rmax; // rmax of the superlayer

    // OUTPUT
    double cell_size;

    double delta_phi;
    int ncells_per_layer; // in one layer
    int ncells_superlayer;

    std::vector<double> start_phi_vec; // start phi [nlayers]
    std::vector<double> signal_wire_r_vec; // signal wire r

    DCSuperLaylerCalculator(int _id, const std::string& _type, int _n, double _rmin, double _rmax)
    : superlayer_id(_id), type(_type), n(_n), rmin(_rmin), rmax(_rmax) {

    }

    void calc() {
        // According to the rmin/rmax and n, calculate the size of cell (inner)
        cell_size = (rmax-rmin) / n;

        // now, let decide the number of cells in each layer.
        // for simplicity, assume they have the same number. 

        double R = rmin;
        delta_phi = cell_size / R;
        ncells_per_layer = 2*M_PI / delta_phi;
        ncells_superlayer = ncells_per_layer*n;

        for (int ilayer = 0; ilayer < n; ++ilayer) {
            double start_phi = 0;
            if (ilayer % 2) {
                start_phi = delta_phi / 2;
            }
            start_phi_vec.push_back(start_phi);

            // signal wire should at center in radial direction
            double r = rmin + 0.5*cell_size + ilayer*cell_size;

            signal_wire_r_vec.push_back(r);
        }

    }

    void dump() {
        std::cout 
        << "SUPER LAYER [" << superlayer_id << "]"
        << " type: " << type
        << " n: " << n
        << " rmin: " << rmin
        << " rmax: " << rmax
        << " cellsize: " << cell_size
        << " ncells_layer: " << ncells_per_layer
        << " ncells_superlayer: " << ncells_superlayer
        << std::endl;

        for (int ilayer = 0; ilayer < n; ++ilayer) {
            std::cout << " - "
            << " ilayer: " << ilayer
            << " r: " << signal_wire_r_vec[ilayer]
            << " startphi : " << start_phi_vec[ilayer]
            << std::endl;
        }
    }

};

static Ref_t create_element(Detector &oddd, xml_h xml,
                                SensitiveDetector sens) {

    xml_det_t x_det = xml;
    std::string detName = x_det.nameStr();

    // Make DetElement
    DetElement driftChamberElement(detName, x_det.id());

    // Make Volume
    // Env, A, U, V
    Volume envVolume;
    Volume wireAVolume;
    Volume wireUVolume;
    Volume wireVVolume;

    if (x_det.hasChild(_U(tubs))) {
        for (xml_coll_t tub(xml, _U(tubs)); tub; ++tub) {
            xml_comp_t x_det_tub = tub;
            std::string shapeName = x_det_tub.nameStr();
            if (shapeName == "DCSupportCylinder") {
                Tube tubeShape(shapeName, x_det_tub.rmin(), x_det_tub.rmax(),
                                x_det_tub.dz());
                Volume tubeVolume(detName, tubeShape,
                                    oddd.material(x_det_tub.materialStr()));
                // tubeVolume.setVisAttributes(oddd, x_det.visStr());
                envVolume = tubeVolume;
            } else if (shapeName == "A") {
                Tube tubeShape(shapeName, x_det_tub.rmin(), x_det_tub.rmax(),
                                x_det_tub.dz());
                Volume tubeVolume(shapeName, tubeShape,
                                    oddd.material(x_det_tub.materialStr()));
                tubeVolume.setVisAttributes(oddd, x_det.visStr());
                wireAVolume = tubeVolume;
            } else if (shapeName == "U") {

            } else if (shapeName == "V") {

            }
        }

    }



    // Place it in the mother
    Volume motherVolume = oddd.pickMotherVolume(driftChamberElement);
    PlacedVolume placedEnv = motherVolume.placeVolume(envVolume);
    placedEnv.addPhysVolID(detName, driftChamberElement.id());
    driftChamberElement.setPlacement(placedEnv);

    // Construct super layers
    std::cout << "START BUILDING SUPER LAYER" << std::endl;
    size_t superlayerNum = 0;
    size_t ncells = 0;
    if (x_det.hasChild(_U(layers))) {
        for (xml_coll_t superlay(xml, _U(layers)); 
            superlay; ++superlay, ++superlayerNum) {
            xml_comp_t x_superlayer = superlay;

            DCSuperLaylerCalculator dcsc(superlayerNum,
                                         x_superlayer.typeStr(),
                                         x_superlayer.number(),
                                         x_superlayer.rmin(),
                                         x_superlayer.rmax());
            dcsc.calc();
            dcsc.dump();

            ncells += dcsc.ncells_superlayer;

            // place the wires
            for (int ilayer = 0; ilayer < dcsc.n; ++ilayer) {
                double start_phi = dcsc.start_phi_vec[ilayer];
                double signal_wire_r = dcsc.signal_wire_r_vec[ilayer];

                for (int icell = 0; icell < dcsc.ncells_per_layer; ++icell) {
                    double phi = start_phi + icell*dcsc.delta_phi;

                    double x = signal_wire_r*cos(phi);
                    double y = signal_wire_r*sin(phi);

                    dd4hep::Transform3D tr(dd4hep::Rotation3D(),
                                           dd4hep::Position(x,y,0));
                    dd4hep::PlacedVolume pv = envVolume.placeVolume(wireAVolume, tr);

                    // debug only
                    if (icell > 10) break;
                }

            }
        }
    }

    std::cout << "Total number of superlayers: " << superlayerNum << std::endl;
    std::cout << "Total number of cells: " << ncells << std::endl;


  return driftChamberElement;

}

DECLARE_DETELEMENT(ODDDriftChamber, create_element)
