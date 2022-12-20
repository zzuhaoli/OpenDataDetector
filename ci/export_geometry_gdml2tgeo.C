void export_geometry_gdml2tgeo(const char* infile, const char* outfile)
{
  // Loading the library and geometry
  gSystem->Load("libGeom");
  TGeoManager::Import(infile);
  gGeoManager->SetVisLevel(4); /// the number here can be changed based on the depth of the visibility level of your gdml file and the detail that you want to visualize it.
  // Export the geometry
  gGeoManager->Export(outfile);
}
