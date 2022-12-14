import argparse
from ROOT import TFile, gROOT, TH1F
import re
from math import pi

TH1F.AddDirectory(False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Scale all TH1 histograms in a file')
    parser.add_argument('--input', '-i', required=True, type=str, nargs='+', help="Names of the input files with results for phi slice")
    parser.add_argument('--output', '-o', required=True, type=str, help="Name of the output file")
    args = parser.parse_args()

    # first determine the axes for phi distributions
    phi_min, phi_max = -180,180
    n_phi_bins = 60
    d_phi = (phi_max - phi_min) / (n_phi_bins -1)
    phi_min -= 0.5 * d_phi
    phi_max += 0.5 * d_phi
    deg2rad = pi/180.

    # now create histograms (phi distributions)
    # get the initial histograms (eta distributions)
    firstfile = TFile(args.input[0],"READ")
    all_hists_eta = {}
    all_hists_phi = {}
    all_hists_phi_normalisation = {}
    for key in firstfile.GetListOfKeys():
        if gROOT.GetClass(key.GetClassName()).InheritsFrom("TH1"):
            h = key.ReadObj()
            all_hists_eta[h.GetName()] = h.Clone()
            all_hists_phi[h.GetName().replace("eta","phi")] = TH1F(h.GetName().replace("eta","phi"), h.GetTitle().replace("eta","phi"), n_phi_bins, phi_min * deg2rad, phi_max * deg2rad)
            all_hists_phi_normalisation[h.GetName().replace("eta","phi")] = TH1F(h.GetName().replace("eta","phi"), h.GetTitle().replace("eta","phi"), n_phi_bins, phi_min * deg2rad, phi_max * deg2rad)
            all_hists_phi[h.GetName().replace("eta","phi")].Sumw2()
            all_hists_phi_normalisation[h.GetName().replace("eta","phi")].Sumw2()

    all_phi_values = []
    phiRegex = re.compile(r'material_budget_phi\-?\d{1,10}\.?\d{0,10}.root')
    for iphi, phi_file in enumerate(args.input):
        phi = int(phiRegex.findall(phi_file)[0][len('material_budget_phi'):-5])
        all_phi_values.append(phi)
        f = TFile(phi_file,"READ")
        for key in f.GetListOfKeys():
            if gROOT.GetClass(key.GetClassName()).InheritsFrom("TH1"):
                h = key.ReadObj()
                # first merge histograms for eta distributions
                if iphi > 0:
                    all_hists_eta[h.GetName()].Add(h)
                # add content to phi distributions
                all_hists_phi[h.GetName().replace("eta","phi")].Fill(all_phi_values[iphi]*deg2rad, h.Integral(1, h.FindBin(1) ) / h.FindBin(1))
                all_hists_phi_normalisation[h.GetName().replace("eta","phi")].Fill(all_phi_values[iphi]*deg2rad)

    # Change names to follow the naming convention of material scan with ACTS
    outfile = TFile(args.output,"RECREATE")
    outfile.cd()
    for hname, h in all_hists_eta.items():
        h.GetXaxis().SetTitle("#eta")
        h.GetYaxis().SetTitle("depth (X_{0})")
        # scale
        h.Scale(1./len(args.input))
        # rename
        if "x0" in hname:
            h.SetName(hname[:-2]+"_x0_vs_eta_all")
        elif "lambda" in hname:
            h.SetName(hname[:-6]+"_l0_vs_eta_all")
        h.Write()

    # Scale phi distributions
    for hkey in all_hists_phi:
        for ibin in range(1, all_hists_phi[hkey].GetXaxis().GetNbins()+1):
            norm = all_hists_phi_normalisation[hkey].GetBinContent(ibin)
            if norm != 0:
                all_hists_phi[hkey].SetBinContent(ibin, all_hists_phi[hkey].GetBinContent(ibin) / norm)
                all_hists_phi[hkey].SetBinError(ibin, all_hists_phi[hkey].GetBinError(ibin) / norm)

    # Change names for phi distributions
    for hname, h in all_hists_phi.items():
        h.SetTitle(h.GetTitle().replace("integrated","integrated for #eta<1"))
        h.GetXaxis().SetTitle("#phi (rad)")
        h.GetYaxis().SetTitle("depth (X_{0})")
        # rename
        if "x0" in hname:
            h.SetName(hname[:-2]+"_x0_vs_phi_all")
        elif "lambda" in hname:
            h.SetName(hname[:-6]+"_l0_vs_phi_all")
        h.Write()

    outfile.Close()
