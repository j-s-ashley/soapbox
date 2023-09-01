/*
mymain05.cc
    AKA Project Tree Rings
Dimuon invariant mass calculation for
multiple events of Z-Boson production, 
decay restricted to muons.
Created 07.20.2023 06:57
Last updated 2023.07.27 06:56
*/

// Headers and Namespaces.
#include "Pythia8/Pythia.h"	// Include Pythia headers.
using namespace Pythia8;	// Let Pythia8:: be implicit.

// Invariant mass calculation
double invMass(double e1, double e2, double x1, double x2, double y1, double y2, double z1, double z2) {
    return sqrt( 
       ( (e1 + e2) * (e1 + e2) ) 
        - ( (x1 + x2) * (x1 + x2) ) 
        - ( (y1 + y2) * (y1 + y2) ) 
        - ( (z1 + z2) * (z1 + z2) ) );
}

// Begin main program.
int main() {

	// Set up generation.
	Pythia pythia;	// Declare a Pythia object.
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on"); // Fermion/anti-fermion collision to photons and Z bosons
    pythia.readString("23:onMode = off"); // Turns off decay for Z ** mayDecay = no for NEVER decay
    pythia.readString("23:onIfAny = -13 13"); // Switch back on Z decays to muons.
    pythia.readString("22:onMode = off"); // Turns off photon decay.
    pythia.readString("22:onIfAny = -13 13"); // Switch back on photon decays to muons.
    pythia.readString("Beams:eCM = 13000."); // 13 TeV CM energy.
    pythia.init(); // Initialize; incoming pp beams is default.

    // Set up histogram
    Hist mZ("N/GeV", 100, 0., 100.); // Name, bins, x-min, x-max.
    HistPlot hpl("mymain05_out"); // File name.
    
    int numEvents = 100000; // Number of events desired.

    // Event loop.
    for(int i = 0; i < numEvents; ++i) {
        if(!pythia.next()) {
            continue; // Error-proofing.
        }
        
        // Particle data.
        int id1 = 0;
        double e1, x1, y1, z1 = 0;
        
        int id2 = 0;
        double e2, x2, y2, z2 = 0;

        // Particle analysis.
        for (int j = 0; j < pythia.event.size(); ++j) {
            // Pick out final particles that are either antimuons or muons.
            if ( (pythia.event[j].isFinal()) && ((pythia.event[j].id() == -13) || (pythia.event[j].id() == 13)) ) {
                if (id1 == 0) { // Checking for first final, skip if found.
                    id1 = pythia.event[j].id();
                    e1 = pythia.event[j].e();
                    x1 = pythia.event[j].px();
                    y1 = pythia.event[j].py();
                    z1 = pythia.event[j].pz();
                } else if (id2 == 0) { // Checking for second final, skip if found.
                    id2 = pythia.event[j].id();
                    e2 = pythia.event[j].e();
                    x2 = pythia.event[j].px();
                    y2 = pythia.event[j].py();
                    z2 = pythia.event[j].pz();
        
                    double m = invMass(e1, e2, x1, x2, y1, y2, z1, z2); // Calculate invariant mass of mother.

                    mZ.fill(m); // Fill histogram.

                    // Reset the variables for the next pair of particles.
                    id1 = e1 = x1 = y1 = z1 = id2 = e2 = x2 = y2 = z2 = 0;
                }
            }
        } // End event record analysis.
    } // Increment event number (i).
    // End event loop.

    // Histogram file name,       title,                                       x-axis label, y-axis label.
    hpl.frame( "DimuonMassDist", "Boson dimuon invariant mass distributions", "mass (GeV)", "events/bin");
    hpl.add( mZ, "h,indigo"); // Type = histogram, color = indigo.
    hpl.plot(true, false); // Plot with a logarithmic y-axis and linear x-axis.

	return 0;	// End main program with error-free return.

}