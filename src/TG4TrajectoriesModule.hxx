#ifndef TG4TrajectoriesModule_hxx_seen
#define TG4TrajectoriesModule_hxx_seen

#include <TAnalysisTruthModuleBase.hxx>
#include <ECaptainSummary.hxx>

#include <TEvent.hxx>
#include <TG4Trajectory.hxx>

#include <TLorentzVector.h>
#include <TClonesArray.h>

#include <set>

namespace CP {
    class TG4TrajectoriesModule;
}


/// Create the G4Trajectories tree in the TruthDir directory to save
/// information about the MC particle paths simulated in the event. The MC
/// information is in a tree called "G4Trajectories" in the "TruthDir"
/// directory of CaptainSummary files.  The tree contains a vector with a
/// TG4Trajectory objects for every particle in the event which passes the
/// module's criteria to be saved.  The tree is defined as:
///
/// * TruthDir/G4Trajectories -- A tree for the true particle trajectories.
///
///   * Trajectory (std::vector<TG4Trajectory>) A branch of TG4Trajectory
///         objects sorted by the trajectory number.
///                 
class CP::TG4TrajectoriesModule : public TAnalysisTruthModuleBase {
public:
    
    class TG4Trajectory;
    
    TG4TrajectoriesModule(const char *name = "G4Trajectories",
                             const char *title = "Summary of G4 trajectories");
    
    virtual ~TG4TrajectoriesModule();
    
    /// This method return kTRUE - the module is always enabled by defult.
    virtual Bool_t IsEnabledByDefault() const {return kTRUE;}
    
    /// Called for the first event, this method checks whether this event is a
    /// real-data event an if so throws an CP::EDataFile(). Otherwise it
    /// returns true.
    virtual Bool_t ProcessFirstEvent(CP::TEvent&);
    
    /// Method for setting behaviour of module. Currently implemented options
    /// are: "saveall" - Save every trajectory longer than the minimum length
    /// regardless of where it occurred
    virtual Bool_t Configure(const std::string &option);
    
    /// Returns the minimum length [mm] of a trajectory before it is saved.
    Double_t GetMinimumTrajectoryLengthToSave() const {return fMinLength;}
		
    /// Set the minimum length [mm] a trajectory must have to be saved by the
    /// module.
    void SetMinimumTrajectoryLengthToSave(Double_t mm) {fMinLength = mm;}
		
    /// Calling this method with kTRUE or no argument specifies that for any
    /// trajectory which would ordinarily be saved, the module should also
    /// save all parent trajectories back to and including the primary
    /// trajectory.
    void SetSaveParentChain(Bool_t yesorno = kTRUE) {
	fSaveParentChain = yesorno;	
    }
		
    /// Initialize Module.  This is where the internal fields of the modules
    /// should be set-up, but don't create the branches.  Modules might not
    /// need to do any thing in this method, but it must be defined.
    virtual void InitializeModule();

    /// Creates the necessary tree and branches for saving the G4
    /// Trajectories information
    virtual void InitializeBranches();
		
    /// Called for each event, this method is the master method for retrieving
    /// and filling the G4 Trajectories information.
    virtual bool FillTree(CP::TEvent&);
		
private:
		
    // Methods
		
    /// Fills a std::set (fSaveList) with the id of every trajectory which
    /// should be saved for the current event.
    void FillSaveList(CP::THandle<CP::TG4TrajectoryContainer> trajectories);
		
    /// Fills the new vectors of entry/exit positions and momenta of the
    /// trajectory for every subdetector the trajectory traverses.
    void FillPoints(CP::TG4Trajectory *const traj,
                    CP::TG4TrajectoriesModule::TG4Trajectory* trajToFill);
		
    /// Returns true if a trajectory needs to be saved, and false oterwise.
    bool SaveTraj(const CP::TG4Trajectory& traj) const;
		
    /// Minimum Length of Trajectories that will be saved in mm.  All primary
    /// particles will be saved regardless of this.
    Double_t fMinLength;

    /// Whether saving a trajectory should also trigger the saving of all the
    /// trajectories in its parent chain.
    Bool_t fSaveParentChain;		

    /// The set of trajectory ids which are to be saved from the current
    /// event.
    std::set<Int_t> fSaveList;		
		
public:
		
    // Members saved by the module to the output tree
		
    /// [TREE MEMBER] std::vector of TG4Trajectory sorted in ascending TrajId
    /// order.
    std::vector<CP::TG4TrajectoriesModule::TG4Trajectory> fTrajectories;

    ClassDef(TG4TrajectoriesModule,1);
};



/// The truth information with an MC particle trajectory in the Trajectory
/// branch of the G4Trajectories tree.  The Trajectory branch of the
/// G4Trajectories tree is a std::vector<TG4Trajectory> object with the
/// trajectories stored in order of increasing TrajId.
class CP::TG4TrajectoriesModule::TG4Trajectory : public TObject {
public:
    
    TG4Trajectory();
    virtual ~TG4Trajectory() {};
		
    /// Trajectory id number.  This uniquely identifies this trajectory within
    /// the event and can be used by other modules to reference trajectories.
    /// When the trajectory was directly generated by a particle from the G4
    /// primary particle stack, then the trajectory id will correspond to the
    /// primary particle track id.  This can lead to some confusion, but the
    /// track id and trajectory id are different concepts (one is for the G4
    /// primary particle, the other is for the G4 trajectory).
    Int_t TrajId;		

    /// Parent particle's TrajId. If this is a trajectory is directly
    /// generated by a primary particle, then the ParentId is zero.  The
    /// TrajId of a trajectory coming directly from a primary particle is the
    /// same as the TrackId of that primary particle. (Be careful, the TrajId
    /// and the TrackId are different concepts, but sometimes get used
    /// interchangable for primary particles.)
    Int_t ParentId;

    /// TrajId of the particle at the top of this particle's parent chain. If
    /// this is a primary trajectory, then the TrackId is equal to TrajId.
    /// Note that for convenience, the particle at the top of the parent chain
    /// may not be a primary particle.  For example, the primary id for a
    /// photon from pizero decay will be the id of the pizero.  Another
    /// example is that the primary id for an electron from a muon decay will
    /// be the id of the parent muon.
    Int_t PrimaryId;

    /// Particle PDG code.  For particles, these are defined by the particle
    /// data group.  For nuclei, these are defined as in PDG2006 with
    /// 10LZZZAAAI [L: total number of strange quarks (usually zero), Z:
    /// charge of nucleus, A: baryon number of nucleus (n+p), I: isomer level
    /// (0 is ground state)]
    Int_t PDG;

    /// Mass of the particle [MeV].
    Float_t Mass;

    /// Charge in units of |e|/3. (e.g. an electron has charge -3)
    Int_t Charge;
		
    /// Vector of TLorentzVector objects storing the positions of the
    /// trajectory points.  The front of the vector is the initial position.
    /// The back of the vector is the final position.
    std::vector<TLorentzVector> Position;

    /// Vector of TVector3 objects that stores the momentum of the of the
    /// particle at each trajectory point.  The front of the vector is the
    /// initial momentum.  The back of the vector is the final momentum.
    std::vector<TVector3> Momentum;

    /// Vector of booleans indicating where the trajectory point is in the
    /// detector.  This is defined by 1nnnppp where vvv is the pre-step
    /// volume, and nnn is the post-step volume.  The pre-step volume is
    /// access with pre = Region[i] % 1000.  The post-step volume is accessed
    /// with post = (Region[i]/1000) % 1000. The sub-detector types are
    /// defined by: 1) the drift region; 100) the cyrostat; 999) Everything
    /// else.  Numbers less than 100 are for active regions.
    std::vector<Int_t> Region;

private:
		
    ClassDef(TG4TrajectoriesModule::TG4Trajectory, 1);
};

#endif
