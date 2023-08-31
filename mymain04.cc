/*
mymain04.cc
    AKA Project Puddle Jumper
Dijet invariant mass calculation for
multiple events of Z-Boson production, 
decay restricted to quark/antiquarks and photons.
Created 07.20.2023 05:36
Last updated 2023.07.28 05:14
*/

// Headers and Namespaces.
#include "Pythia8/Pythia.h"	// Include Pythia headers.
#include "fastjet/ClusterSequence.hh" // Include fastjet headers.
using namespace fastjet;	// Let fastjet:: be implicit.    
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
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off"); // Turns off decay for Z ** mayDecay = no for NEVER decay
    pythia.readString("23:onIfAny = 1 2 3 4 5"); // Switch back on Z decays to quarks
    pythia.readString("22:onMode = off"); // Turns off photon decay.
    pythia.readString("22:onIfAny = 1 2 3 4 5"); // Switch back on photon decays to quarks.
    pythia.readString("Beams:eCM = 13000."); // 13 TeV CM energy.
    pythia.init(); // Initialize; incoming pp beams is default.

    // Set up histogram
    Hist pTZ("N/GeV", 100, 0., 100.); // Data description, bins, x-min, x-max.
    HistPlot hpl("mymain04_out"); // File name.
    
    int numEvents = 100000; // Number of events desired.
    vector<PseudoJet> particles; // FastJet vector for particle 4-momentum.
    
    double R = 0.4; // Set jet radius.
    JetDefinition jet_def(antikt_algorithm, R); // Define jet algorithm.

    // Event and jet loop.
    for(int i = 0; i < numEvents; ++i) {
        if(!pythia.next()) {
            continue; // Error-proofing.
        }
        
        for(int j = 0; j < pythia.event.size(); ++j) {      // Analyze event record looking for
            if(pythia.event[j].isFinal()) {                 // final particles and push their 
                  particles.push_back( PseudoJet(            // 4-momentum to the particles vector.
                      pythia.event[j].px(), pythia.event[j].py(),
                      pythia.event[j].pz(), pythia.event[j].e()) 
                  );
            } // End event record analysis.
        }
        
        // Cluster and extract jets
        ClusterSequence cs(particles, jet_def);
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); // Jets is overwritten each iteration.
        
        // Send data to invariant mass function
        double m = invMass( 
            jets[0].e(), jets[1].e(),
            jets[0].px(), jets[1].px(),
            jets[0].py(), jets[1].py(),
            jets[0].pz(), jets[1].pz()
        );
        pTZ.fill(m); // Log invariant mass.

        particles.clear(); // Empty particles vector.
    } // Increment event number (i).
    // End event and jet loop.
    
    hpl.frame( "DijetMassDist", "Boson Invariant Mass Distributions", "m (GeV)", "events/bin"); 
    hpl.add( pTZ, "h,indigo");
    hpl.plot(true, false);

	return 0;	// End main program with error-free return.

}