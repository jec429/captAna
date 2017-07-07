#ifndef TG4VerticesModule_hxx_seen
#define TG4VerticesModule_hxx_seen

#include <TAnalysisTruthModuleBase.hxx>
#include <ECaptainSummary.hxx>

#include <TEvent.hxx>
#include <TG4Trajectory.hxx>

#include <TLorentzVector.h>
#include <TClonesArray.h>

#include <set>

namespace CP {
    class TG4VerticesModule;
}


/// Create the G4Vertices tree in the TruthDir directory to save information
/// about the MC primary vertices and particles.  The MC informaion is in a
/// tree called "G4Vertices" in the "TruthDir" directory of the DST file. The
/// tree contains the Vertex vector with a TG4Vertex object for every vertex,
/// and the Particle vector of TG4Particle for every particle attached to a
/// vertex.  The tree is defined as:
///
/// * TruthDir/G4Vertices -- A tree for the primary particles started in G4.
///
///   * Vertex (std::vector<TG4Vertex>) A branch of TG4Vertex objects sorted
///         by the vertex number.
///
///   * Particle (std::vector<TG4Particle>) A branch of TG4Particle objects
///         sorted by primary particle id.
///                 
class CP::TG4VerticesModule : public TAnalysisTruthModuleBase {
public:
    
    class TG4Vertex;
    class TG4Particle;
    
    TG4VerticesModule(const char *name = "G4Vertices",
                             const char *title = "Summary of G4 trajectories");
    
    virtual ~TG4VerticesModule();
    
    /// Called for the first event, this method checks whether this event is a
    /// real-data event an if so throws an CP::EDataFile(). Otherwise it
    /// returns true.
    virtual Bool_t ProcessFirstEvent(CP::TEvent&);
    
    /// Initialize Module.  This is where the internal fields of the modules
    /// should be set-up, but don't create the branches.  Modules might not
    /// need to do any thing in this method, but it must be defined.
    virtual void InitializeModule();

    /// Creates the necessary tree and branches for saving the Truth
    /// Vertices information
    virtual void InitializeBranches();
		
    /// Called for each event, this method is the master method for retrieving
    /// and filling the Truth Vertices information.
    virtual bool FillTree(CP::TEvent&);

    /// [TREE MEMBER] std::vector of TG4Vertex for each vertex in the event.
    std::vector<CP::TG4VerticesModule::TG4Vertex> GetVerticies() {
      return fVertices;
    }

    std::vector<CP::TG4VerticesModule::TG4Particle> GetParticles() {
      return fParticles;
    }
    
		
private:
		
    /// [TREE MEMBER] std::vector of TG4Vertex for each vertex in the event.
    std::vector<CP::TG4VerticesModule::TG4Vertex> fVertices;

    /// [TREE MEMBER] std::vector of TG4Particle for each particle attached to a
    /// vertex in the event.
    std::vector<CP::TG4VerticesModule::TG4Particle> fParticles;

    ClassDef(TG4VerticesModule,1);
};

/// The Vertex branch of the G4Vertices tree contains the truth information
/// associated with a vertex from the Monte Carlo simulations.
class CP::TG4VerticesModule::TG4Vertex : public TObject {
public:
    
    TG4Vertex();
    virtual ~TG4Vertex() {};

    /// Vertex id number.  This uniquely identifies this vertex within
    /// the event and can be used by other modules to reference a vertex.
    Int_t VertexId;

    /// The vertex position.
    TLorentzVector Position;

    /// A vector of id numbers for the truth particles associated with this
    /// vertex.
    std::vector<Int_t> Particles;

private:
		
    ClassDef(TG4VerticesModule::TG4Vertex, 1);
};

/// The Particle branch of the G4Vertices tree contains the truth information
/// associated with a particle from the Monte Carlo simulations.  
class CP::TG4VerticesModule::TG4Particle : public TObject {
public:
    
    TG4Particle();
    virtual ~TG4Particle() {};

    /// Track id number.  This uniquely identifies this particle within the
    /// event and can be used by other modules to reference the particle.  The
    /// track id also corresponds to the trajectory that was created by this
    /// particle.
    Int_t TrackId;

    // add parent
    Int_t ParentId;

    /// The vertex id that this particle comes from.
    Int_t VertexId;

    /// Particle PDG code.  For particles, these are defined by the particle
    /// data group.  For nuclei, these are defined as in PDG2006 with
    /// 10LZZZAAAI [L: total number of strange quarks (usually zero), Z:
    /// charge of nucleus, A: baryon number of nucleus (n+p), I: isomer level
    /// (0 is ground state)]
    Int_t PDG;

    /// The particle mass
    Float_t Mass;
    
    /// The particle charge in units of |e|/3, so an electron has a charge of
    /// -3.
    Int_t Charge;

    /// The momentum of the particle
    TLorentzVector Momentum;

private:
		
    ClassDef(TG4VerticesModule::TG4Particle, 1);
};

#endif
