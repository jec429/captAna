//c++
#include <iostream>
#include <vector>
#include <assert.h>
#include <memory>
#include <numeric>
#include <cmath>
#include <complex>
// captain
#include "dstUtils.hxx"
#include <HEPUnits.hxx>
#include <HEPConstants.hxx>
#include <TRootOutput.hxx>
#include <TDigitContainer.hxx>
#include <TCalibPulseDigit.hxx>
#include <TPulseDigit.hxx>
#include <TPulseHit.hxx>
#include <TReconHit.hxx>
#include <TChannelId.hxx>
#include <TChannelCalib.hxx>
#include <TG4PrimaryVertex.hxx>
#include <TG4PrimaryParticle.hxx>
#include <TMCChannelId.hxx>
#include <TMCDigit.hxx>

//
#include <TEvent.hxx>
#include <TEventContext.hxx>
#include <CaptGeomId.hxx>
#include <TChannelInfo.hxx>
#include <TGeometryInfo.hxx>
#include "CaptGeomIdDef.hxx"
#include <TGeomIdManager.hxx>
#include <TManager.hxx>
#include <TReconTrack.hxx>
#include <TTrackState.hxx>
#include <TReconCluster.hxx>
#include "TClusterState.hxx"
#include <TReconShower.hxx>
#include <TElectronicsResponse.hxx>
#include <TTPCChannelId.hxx>
#include "TG4TrajectoriesModule.hxx"
#include "TG4VerticesModule.hxx"
#include <TG4Trajectory.hxx>
#include <TG4HitSegment.hxx>
#include <HEPUnits.hxx>
#include <TG4HitSegment.hxx>
#include <TUnitsTable.hxx>
#include <HEPUnits.hxx>
#include <eventLoop.hxx>



// root

#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>


#include <TVirtualFFT.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TFitResult.h>
#include <TTree.h>
// captAna
#include <TCapEvent.hxx>
#include <TCapTruth.hxx>
#include <TCapTrack.hxx>


namespace CP {
    class captAna;
}

using namespace std;

class CP::captAna: public CP::TEventLoopFunction {
  public:
    captAna() {}
    virtual ~captAna() {};
    void Initialize(void);
    bool operator () (CP::TEvent& event);
    void Finalize(CP::TRootOutput * const output);
    void InitTruth();
    typedef TMatrixTSym<float> MomentMatrix;
    

  private:
    std::string ChIdString( CP::TTPCChannelId tpcid ) {
      char chname[80]; 
      sprintf(chname,"tpcid_%u_%u_%u",tpcid.GetCrate(),tpcid.GetFEM(),tpcid.GetChannel());
      return std::string(chname);
    }

    CP::TDigit* GetDigit( CP::THandle<CP::TDigitContainer> container , CP::TChannelId chID ){

      if(!container) {
        printf(" container pointer  is null \n");
        return NULL;
      }
      //printf(" searching in  container of size %i \n",container->size());

      for (CP::TDigitContainer::iterator td = container->begin(); td != container->end(); ++td) {
        //cout << " looking for  " << chID.AsString() << " ?? " << (*td)->GetChannelId().AsString() << endl;
        if( (*td)->GetChannelId()==chID) return (*td);
      }
      return NULL;

    }

    // find the associated cluster hit
    int findCluster(CP::THandle<CP::THit> phit) { 
      CP::TChannelId chid = phit->GetChannelId();
      CP::TGeometryId geoid = CP::TChannelInfo::Get().GetGeometry(chid);
      int plane = CP::GeomId::Captain::GetWirePlane(geoid);
      int wire = (int) CP::GeomId::Captain::GetWireNumber(geoid) + 0.5;
      CP::TTPCChannelId tpcid(chid);
      int nChannel =tpcid.GetChannel();
      int fem = tpcid.GetFEM();
      //printf(" looking for %i %i %i %i \n", nChannel, fem, plane, wire );
      int icluster = -1;
      for( unsigned ihit =0; ihit < tcapevent->hits.size() ; ++ ihit) {
        CP::TCapHit chit = tcapevent->hits[ihit];
        if( chit.nchan == nChannel && chit.fem == fem && chit.plane == plane &&  chit.wire == wire)
          icluster = ihit;
      }
      return icluster;
    }



    //std::string fLSOption;
    //TH1D* fWires[2000];
    bool isFirstEvent;
    enum {MAXHIST=2};
    enum {MAXFEM=16};
    enum {MAXASIC=16};
    enum {NSIZE=9596};
    CP::TElectronicsResponse* fElectronicsResponse;
    int nWaveHist;
    CP::TTPCChannelId tpcid_list[MAXHIST];
    int event_list[MAXHIST];
    TTree* tdrift;
    TH1D* hWaveFit;
    TH1D* hWaveHits[MAXHIST];
    TH1D* hWaveSamples[MAXHIST];
    TH1D* hDWaveSamples[MAXHIST];
    TH1D* hFFTSamples[MAXHIST];
    TH1D* hElecResponse[MAXHIST];
    TH1D* hFFT_FEM[MAXFEM][MAXASIC]; 
    TH1D* hFFT_DEC[MAXFEM][MAXASIC]; 
    TNtuple *ntPmt;
    // cluster histos
    TH1D* hClusterQ;
    int maxadc;
    unsigned total_events;
    unsigned ilate;
    //
    typedef struct{
      int iev;
      int crate;
      int fem;
      int asic;
      int nchan;
      int plane;
      int wire;
      float baseMedian;
      float baseSigma;
      float sampleSigma;
      float fitMean;
      float fitMeanErr;
      float fitSlope;
      float fitSlopeErr;
      float fitChi; // per dof 
      float npower;
    } WAVE;
    WAVE rawWave;
  
     //bool fQuiet;
    //bool fPrintList;

    /// The fft class to take the fourier transform.
    TVirtualFFT *fFFT;
    TVirtualFFT *fFFTResponse;

    /// The fft class to take the inverse fourier transform.
    TVirtualFFT *fInverseFFT;

    TTree *anaTree;
    CP::TCapEvent *tcapevent;
    CP::TCapTruth *tcaptruth;
    CP::TCapCluster *tcalibcluster;

    // from captSummary
    CP::TG4TrajectoriesModule* fTG4TrajectoriesModule;
    CP::TG4VerticesModule* fTG4VerticiesModule;

    // 
    vector<std::string> reconClassName;
    typedef enum {RTRACK=0,RSHOWER,RCLUSTER,RVERTEX,RPID,UNKNOWN} reconType;
    int typeList[UNKNOWN+1];

    // cache hits temporarily
    std::vector<CP::TCapHit> fhits; 
    std::vector<CP::TCapHitSegment> fsegments; 

    // fill methods
    void fillCapHit(CP::TCapHit* capHit,CP::THandle<CP::TPulseHit> hit, double noise);
    void fillCapReconHit(CP::TCapReconHit* capHit,CP::THandle<CP::TReconHit> hit);
    void fillCapTrack(CP::TCapTrack* track, CP::THandle<CP::TReconTrack> rtrack);
    void fillCapCluster(CP::TCapCluster* cluster, CP::THandle<CP::TReconCluster> rcluster);    
    void fillCapShower(CP::TCapShower* shower, CP::THandle<CP::TReconShower> rshower);    
    void fillTrackNode(CP::TCapTrackNode* node, CP::THandle<CP::TTrackState> tstate);
    void fillTruth(CP::TEvent& event);
    void fillTruthPoints(CP::TCapTrajectory *atraj, CP::TG4Trajectory* traj, int& ndrift, int& ncryo, int &nother );// return number in drift region

    void fillTruthParticle(CP::TCapPrimaryParticle *aparticle, CP::TG4PrimaryParticle* p , int ivertex); 
     

    void printTypeList() {
      printf(" ************  list of types  *****  \n");
      printf(" track   %i \n ",typeList[RTRACK]);
      printf(" shower  %i \n ",typeList[RSHOWER]);
      printf(" cluster %i \n ",typeList[RCLUSTER]);
      printf(" vertex  %i \n ",typeList[RVERTEX]);
      printf(" pid     %i \n ",typeList[RPID]);
      printf(" unknown %i \n ",typeList[UNKNOWN]);
    }
    
    
};

