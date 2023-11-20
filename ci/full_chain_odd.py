#!/usr/bin/env python3
import argparse
import pathlib, acts, acts.examples
import os, argparse, pathlib, acts, acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    addCKFTracks,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addAmbiguityResolutionML,
    AmbiguityResolutionMLConfig,
    addVertexFitting,
    VertexFinder,
)
from acts.examples.odd import getOpenDataDetector

def getOpenDataDetector(odd_dir, mdecorator=None):
    import acts.examples.dd4hep
    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[str(odd_dir / "xml/OpenDataDetector.xml")]
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()
    if mdecorator is None:
        mdecorator = acts.JsonMaterialDecorator(
            rConfig=config,
            jFileName=str(odd_dir / "config/odd-material-mapping-config.json"),
            level=acts.logging.WARNING,
        )

    trackingGeometry, deco = detector.finalize(dd4hepConfig, mdecorator)

    return detector, trackingGeometry, deco

parser = argparse.ArgumentParser(description="OpenDataDetector full chain example")
parser.add_argument(
    "-n", "--events", type=int, default=100, help="Number of events to run"
)
parser.add_argument(
    "-s", "--skip", type=int, default=0, help="Number of events to skip"
)
parser.add_argument(
    "-j",
    "--jobs",
    type=int,
    default=-1,
    help="Number of threads to use. Default: -1 i.e. number of cores",
)
parser.add_argument(
    "-o",
    "--output",
    type=pathlib.Path,
    default=pathlib.Path.cwd(),
    help="Output directories. Default: $PWD",
)

args = vars(parser.parse_args())

u = acts.UnitConstants
geoDir = pathlib.Path(__file__).parent.parent
outputDir = args["output"]
# acts.examples.dump_args_calls(locals())  # show python binding calls

oddMaterialMap = geoDir / "data/odd-material-maps.root"
oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = geoDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector, trackingGeometry, decorators = getOpenDataDetector(
    geoDir, mdecorator=oddMaterialDeco
)
field = detector.field
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args["events"],
    numThreads=args["jobs"],
    outputDir=str(outputDir),
    trackFpes=False,
)

addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-4.0, 4.0, True),
        ParticleConfig(2, acts.PdgParticle.eMuon, True),
        rnd=rnd,
)

addFatras(
    s,
    trackingGeometry,
    field,
    preSelectParticles=ParticleSelectorConfig(),
    enableInteractions=True,
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
    rnd=rnd,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
    rnd=rnd,
)

addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(),
    geoSelectionConfigFile=oddSeedingSel,
    outputDirRoot=outputDir,
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(
        pt=(0.0, None),
        absEta=(None, 3.0),
        loc0=(-4.0 * u.mm, 4.0 * u.mm),
        nMeasurementsMin=7,
    ),
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
)

addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(
        maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
    ),
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
)

addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.Iterative,
    outputDirRoot=outputDir,
)

s.run()
