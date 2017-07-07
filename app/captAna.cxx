#include <eventLoop.hxx>
#include <captAna.hxx>


 
 void CP::captAna::fillTrackNode(CP::TCapTrackNode* node,CP::THandle<CP::TTrackState> tstate) {
  if(!tstate) return;
  node->edep= tstate->GetEDeposit();
  node->edepVariance= tstate->GetEDepositVariance();
  node->positionVec = tstate->GetPosition();
  node->positionVariance = tstate->GetPositionVariance();
  return ;
}

void CP::captAna::fillCapHit(CP::TCapHit* capHit, CP::THandle<CP::TPulseHit> hit, double noise=0) {
  CP::TChannelId chid = hit->GetChannelId();
  CP::TGeometryId geoid = CP::TChannelInfo::Get().GetGeometry(chid);
  int plane = CP::GeomId::Captain::GetWirePlane(geoid);
  double wire = CP::GeomId::Captain::GetWireNumber(geoid) + 0.5;
  CP::TTPCChannelId tpcid(chid);
  int nChannel =tpcid.GetChannel();
  int nASIC= nChannel>>4; 
  capHit->nchan= nChannel;
  capHit->fem= tpcid.GetFEM();
  capHit->crate=tpcid.GetCrate();
  capHit->asic=nASIC;
  capHit->plane=plane;
  capHit->wire=int(wire);
  capHit->q = hit->GetCharge();
  capHit->dq = hit->GetChargeUncertainty();
  capHit->x = hit->GetPosition().X();
  capHit->dx = hit->GetUncertainty().X();
  capHit->y = hit->GetPosition().Y();
  capHit->dy = hit->GetUncertainty().Y();
  capHit->z = hit->GetPosition().Z();
  capHit->dz = hit->GetUncertainty().Z();
  capHit->t = hit->GetTime();
  capHit->dt = hit->GetTimeUncertainty();
  capHit->rmst = hit->GetTimeRMS();
  capHit->noise = noise;
  capHit->tstart = hit->GetTimeStart();
  capHit->tstop = hit->GetTimeStop();
  int nsamples = hit->GetTimeSamples();
  float qpeak=0.0;
  for(int isample =0; isample < nsamples; ++isample) {
    float qsample = (float) hit->GetTimeSample(isample);
    if(qsample>qpeak) qpeak=qsample;
    capHit->qsamples.push_back(qsample);  
  }
  capHit->qpeak=qpeak;
  capHit->nsamples = int(capHit->qsamples.size());
  //capHit->print(int(wire));
}

// recon hit

void CP::captAna::fillCapReconHit(CP::TCapReconHit* capHit, CP::THandle<CP::TReconHit> hit) {
  CP::TChannelId chid = hit->GetChannelId();
  CP::TGeometryId geoid = CP::TChannelInfo::Get().GetGeometry(chid);
  int plane = CP::GeomId::Captain::GetWirePlane(geoid);
  double wire = CP::GeomId::Captain::GetWireNumber(geoid) + 0.5;
  CP::TTPCChannelId tpcid(chid);
  int nChannel =tpcid.GetChannel();
  int nASIC= nChannel>>4; 
  capHit->nchan= nChannel;
  capHit->fem= tpcid.GetFEM();
  capHit->crate=tpcid.GetCrate();
  capHit->asic=nASIC;
  capHit->plane=plane;
  capHit->wire=int(wire);
  capHit->q = hit->GetCharge();
  capHit->dq = hit->GetChargeUncertainty();
  capHit->x = hit->GetPosition().X();
  capHit->dx = hit->GetUncertainty().X();
  capHit->y = hit->GetPosition().Y();
  capHit->dy = hit->GetUncertainty().Y();
  capHit->z = hit->GetPosition().Z();
  capHit->dz = hit->GetUncertainty().Z();
  capHit->t = hit->GetTime();
  capHit->dt = hit->GetTimeUncertainty();
  capHit->rmst = hit->GetTimeRMS();
  capHit->position = hit->GetPosition();
  capHit->positionUnc= hit->GetUncertainty();
  capHit->positionRMS= hit->GetRMS();
  capHit->nconstituents= hit->GetConstituentCount();
  for(int icon =0; icon < capHit->nconstituents; ++icon){
    CP::THandle<CP::THit> hcon = hit->GetConstituent(icon);
    int icluster = findCluster(hcon);
    if(icluster>=0) {
      capHit->conhits.push_back(icluster);
    } else {
      printf(" warning-- constintuent hit not found \n");
      hcon->ls();
    }
  }
    
  capHit->ndigits= hit->GetDigitCount();
}




//track
void CP::captAna::fillCapTrack(CP::TCapTrack* track, CP::THandle<CP::TReconTrack>  rtrack) {
  track->edep = rtrack->GetEDeposit();
  track->positionVec= rtrack->GetPosition();
  track->positionVariance= rtrack->GetPositionVariance();
  track->dimensions = rtrack->GetDimensions();
  track->isXTrack = rtrack->IsXTrack();
  track->isYTrack = rtrack->IsYTrack();
  track->isZTrack = rtrack->IsZTrack();
  track->directionVec = rtrack->GetDirection();
  track->mass = rtrack->GetMass();
  track->width = rtrack->GetWidth();
  track->quality = rtrack->GetQuality();
  track->ndof = rtrack->GetNDOF();
// from state, front and back
   
  // fill hits for this track
  CP::THandle<CP::THitSelection> rhits = rtrack->GetHits();
  for (CP::THitSelection::iterator rhit = rhits->begin(); rhit != rhits->end(); ++rhit) {
    CP::TCapReconHit capHit;
    CP::THandle<CP::TReconHit> theHit = *rhit;
    if(theHit) {
      fillCapReconHit(&capHit,theHit);
      track->trkhits.push_back(capHit); 
    } 
  }
  // fill nodes 
  CP::TReconNodeContainer nodeContainer = rtrack->GetNodes();  
  for(TReconNodeContainer::iterator nitr =nodeContainer.begin(); nitr !=nodeContainer.end(); ++nitr ) {
    CP::THandle<CP::TTrackState> tstate = (*nitr)->GetState();
    if(tstate) {
      CP::TCapTrackNode node;
      fillTrackNode(&node,tstate);
      track->nodes.push_back(node);
    }
  }
  /* 
    CP::THandle<CP::TTrackState> frontState = rtrack->GetState();
    if(frontState) {
      track->frontPos   = frontState->GetPosition();
      track->frontVar   = frontState->GetPositionVariance();
      track->frontDir   = frontState->GetDirection().Unit();
      track->frontDvar   = frontState->GetDirectionVariance();
    }
    CP::THandle<CP::TTrackState> backState = rtrack->GetState();
    if(backState) {
      track->backPos   = backState->GetPosition();
      track->backVar   = backState->GetPositionVariance();
      track->backDir   = backState->GetDirection().Unit();
      track->backDvar   = backState->GetDirectionVariance();
    }
    */

    TLorentzVector rdiff = track->nodes.front().positionVec - track->nodes.back().positionVec;
    track->length = sqrt(rdiff.X()*rdiff.X()+rdiff.Y()*rdiff.Y()+rdiff.Z()*rdiff.Z());

  
}
//cluster
void CP::captAna::fillCapCluster(CP::TCapCluster* cluster, CP::THandle<CP::TReconCluster> rcluster) {
   // gotta put this in CP::THandle<CP::TClusterState> GetState() const {
    cluster->edep= rcluster->GetEDeposit();
    cluster->edepVariance=rcluster->GetEDepositVariance();
    cluster->positionVec=rcluster->GetPosition();
    cluster->positionVariance=rcluster->GetPositionVariance();
    cluster->dimensions = rcluster->GetDimensions();
    cluster->isXCluster = rcluster->IsXCluster();
    cluster->isYCluster = rcluster->IsYCluster();
    cluster->isZCluster = rcluster->IsZCluster();
    MomentMatrix mmatrix = rcluster->GetMoments();
    // why do these not link?
    //cluster->longAxis =rcluster->GetLongAxis();
    //cluster->majorAxis =rcluster->GetMajorAxis();
    //cluster->minorAxis =rcluster->GetMinorAxis();
    //cluster->longExtent =rcluster->GetLongExtent();
    //cluster->majorExtent =rcluster->GetMajorExtent();
    //cluster->minorExtent =rcluster->GetMinorExtent();

    CP::TReconCluster::MomentMatrix rmoments = rcluster->GetMoments();
    //cluster->fMoments = rmoments;
    for(int im =0; im<3 ; ++im) {
      for (int jm=0; jm<3; ++jm ) cluster->moments[im][jm] = rmoments(im,jm);
    }
    // add hits for the cluster 
    CP::THandle<CP::THitSelection> rhits = rcluster->GetHits();
    for (CP::THitSelection::iterator rhit = rhits->begin(); rhit != rhits->end(); ++rhit) {
      CP::TCapReconHit capHit;
      CP::THandle<CP::TReconHit> theHit = *rhit;
      if(theHit) {
        fillCapReconHit(&capHit,theHit);
        cluster->hits.push_back(capHit); 
      }
    }
 } 

// shower
void CP::captAna::fillCapShower(CP::TCapShower* shower, CP::THandle<CP::TReconShower> rshower) {
  // gotta put this in CP::THandle<CP::TShowerState> GetState() const {
  //printf(" filling shower \n");
  shower->edep= rshower->GetEDeposit();
  shower->positionVec=rshower->GetPosition();
  shower->positionVariance=rshower->GetPositionVariance();
  shower->dimensions = rshower->GetDimensions();
  shower->isXShower = rshower->IsXShower();
  shower->isYShower = rshower->IsYShower();
  shower->isZShower = rshower->IsZShower();
  shower->coneAngle = rshower->GetConeAngle();
  //printf(" %f %f \n",shower->edep,shower->coneAngle);
  // add hits for the shower 
  CP::THandle<CP::THitSelection> rhits = rshower->GetHits();

  for (CP::THitSelection::iterator rhit = rhits->begin(); rhit != rhits->end(); ++rhit) {
    CP::TCapReconHit capHit;
    CP::THandle<CP::TReconHit> theHit = *rhit;
    if(theHit) {
      fillCapReconHit(&capHit,theHit);
      shower->swrhits.push_back(capHit); 
    }
  }
}

void CP::captAna::InitTruth(){
  // truth branch in anaTree
  tcaptruth= new CP::TCapTruth();
  anaTree->Branch("truth",&tcaptruth);
  
  gDirectory->cd("/");
  fTG4TrajectoriesModule= new CP::TG4TrajectoriesModule(); 
  fTG4VerticiesModule = new CP::TG4VerticesModule();      
  std::string dirname = fTG4TrajectoriesModule->GetDirectoryName();
  gDirectory->cd("/");
  TObject *dirExist;
  dirExist = gDirectory->Get(dirname.c_str());
  if(!dirExist) gDirectory->mkdir(dirname.c_str());
  gDirectory->cd(dirname.c_str());
  fTG4TrajectoriesModule->Initialize(new TTree(fTG4TrajectoriesModule->GetName(), fTG4TrajectoriesModule->GetTitle()));


  dirname = fTG4VerticiesModule->GetDirectoryName();
  gDirectory->cd("/");
  dirExist = gDirectory->Get(dirname.c_str());
  if(!dirExist) gDirectory->mkdir(dirname.c_str());
  gDirectory->cd(dirname.c_str());
  fTG4VerticiesModule->Initialize(new TTree(fTG4VerticiesModule->GetName(), fTG4VerticiesModule->GetTitle()));
}

void CP::captAna::fillTruthParticle(CP::TCapPrimaryParticle *aparticle, CP::TG4PrimaryParticle* p , int ivertex ) {
  aparticle->TrackId = p->GetTrackId();
  aparticle->ParentId = p->GetParentId();
  aparticle->VertexId = ivertex;
  aparticle->PDG = p->GetPDGCode();
  //printf(" vertex %i track %i parent %i (%i) pdg %i \n",aparticle->VertexId,aparticle->TrackId,p->GetParentId(),aparticle->ParentId,aparticle->PDG); 
  TParticlePDG* pdgPart = TDatabasePDG::Instance()->GetParticle(aparticle->PDG);
  // make sure particle exists in data base
  if(pdgPart) {
    aparticle->Mass = CP::utils::FindPDGMass(aparticle->PDG);
    aparticle->Charge = CP::utils::FindPDGCharge(aparticle->PDG);
  }
  aparticle->Momentum = p->GetMomentum();
  //printf(" blah \n");
  //p->ls();
  //aparticle->print(p->GetTrackId());
}

// from TG4TrajectoriesModule  
// I put this in because I wanted to be sure the correctly ordered trajectory was filled
void CP::captAna::fillTruthPoints(CP::TCapTrajectory *atraj, CP::TG4Trajectory* traj, int& ndrift, int& ncryo, int& nother) {
    ndrift=0;
    ncryo=0;
    nother=0;
    // Get the points from the trajectory.
    const CP::TG4Trajectory::Points& points = traj->GetTrajectoryPoints();
    if (points.empty()) {
        printf("Found TG4Trajectory with no trajectory points");
        return ;
    }
    for (CP::TG4Trajectory::Points::const_iterator p = points.begin(); p != points.end(); ++p) {
        std::string volumeName = p->GetVolumeName();
        // Figure out the volume number for this.  Be careful with the
        // conditions since the first one that's met determine the volume number.
        int thisVolume = 0;
        if (volumeName.find("Drift") != std::string::npos) {
            thisVolume = 1;
        }
        else if (volumeName.find("Cryostat") != std::string::npos) {
            thisVolume = 50;
        }
        else {
            thisVolume = 999;
        }

        // Find the next volume name.
        CP::TG4Trajectory::Points::const_iterator n = p;
        int nextVolume = 0;
        if (++n != points.end()) {
            volumeName = n->GetVolumeName();
            if (volumeName.find("Drift") != std::string::npos) {
                nextVolume = 1;
            }
            else if (volumeName.find("Cryostat") != std::string::npos) {
                nextVolume = 50;
            }
            else {
                nextVolume = 999;
            }
        }

        atraj->Position.push_back(p->GetPosition());
        atraj->Momentum.push_back(p->GetMomentum());
        atraj->Region.push_back(1000000 + 1000*nextVolume + thisVolume);
        if(thisVolume==1) ++ndrift;
        else if(thisVolume==50) ++ncryo;
        else  ++nother;
    }    
}


void CP::captAna::fillTruth(CP::TEvent& event)
{
  // already done in processEvent setep
  //fTG4VerticiesModule->FillTree(event); 
  //fTG4VerticiesModule->GetOutputTree()->Fill(); 
  //fTG4TrajectoriesModule->FillTree(event);  
  //fTG4TrajectoriesModule->GetOutputTree()->Fill();

  tcaptruth->RunId=event.GetRunId();
  tcaptruth->EventId=event.GetEventId();
  //tcapTruth->SubrunId=event->GetSubrunId();

  //add additional verticies
  CP::THandle<CP::TG4PrimaryVertexContainer> vertices
    =   event.Get<CP::TG4PrimaryVertexContainer>("truth/G4PrimVertex00");
  
  //printf(" number of verticies in TG4PrimaryVertexContainer %i \n",(int) vertices->size());

  for (CP::TG4PrimaryVertexContainer::iterator v = vertices->begin();  v != vertices->end(); ++v) {
    CP::TCapPrimaryVertex avert;
    avert.VertexId =  (int) tcaptruth->pvertex.size();
    avert.Position = v->GetPosition();
    avert.Generator = TString(v->GetGeneratorName());
    avert.Reaction = TString(v->GetReaction());
    printf(" vertex id %i \n",avert.VertexId);

    // add starting particles
    TG4PrimaryVertexContainer infoVertex =  v->GetInfoVertex();
    if(infoVertex.size()<1) {
      printf(" InfoVertex size is 0 \n"); 
    } else {
      CP::TG4PrimaryParticleContainer infoParticles = infoVertex[0].GetPrimaryParticles();
      printf(" infopartiles  %i \n\n\n", (int) infoParticles.size());
      avert.Ninfo = (int)  infoParticles.size();
      for (unsigned infop =0 ; infop < infoParticles.size() ; ++infop) {
        CP::TCapPrimaryParticle aparticle;
        printf(" info particle  %i\n",(int) infop);
        infoParticles[infop].ls();
        fillTruthParticle(&aparticle, &infoParticles[infop], avert.VertexId);
        avert.Particles.push_back(tcaptruth->pparticle.size());   
        tcaptruth->pparticle.push_back(aparticle);  
        //printf(" size of pparticle %i \n",(int) tcaptruth->pparticle.size()); 
      } 
    }
   
    //printf(" infovertex %i \n\n\n", (int) infoVertex.size());
    //printf(" gen %s reaction %s \n",infoVertex[0].GetGeneratorName(), infoVertex[0].GetReaction());
               
    //fill in rest of particles 
    CP::TG4PrimaryParticleContainer& particles = v->GetPrimaryParticles();
    for (unsigned ipart =0 ; ipart < particles.size() ; ++ipart) {
      CP::TCapPrimaryParticle aparticle;
      fillTruthParticle(&aparticle, &particles[ipart], avert.VertexId );
      avert.Particles.push_back(tcaptruth->pparticle.size());   
      tcaptruth->pparticle.push_back(aparticle);             
    }
    tcaptruth->pvertex.push_back(avert);
    printf(" size of pvertex %i  pparticle %i \n",(int) tcaptruth->pvertex.size(),(int) tcaptruth->pparticle.size());        
  }


  // now truth hits
  CP::THandle<CP::TDataVector> truthHits = event.Get<CP::TDataVector>("truth/g4Hits");
  CP::THandle<CP::TG4TrajectoryContainer> fTrajectories = event.Get<CP::TG4TrajectoryContainer>("truth/G4Trajectories");
  
  fsegments.clear();
  if(!truthHits) {
    printf("truth/g4Hits not found !!! ");
  } else {
    //printf("xxxxxxx truth/g4Hits found xxxxx \n ");
    
    //double minEnergy = 0.18*unit::MeV/unit::mm;
    //double maxEnergy = 3.0*unit::MeV/unit::mm;
    for (CP::TDataVector::iterator h = truthHits->begin(); h != truthHits->end(); ++h) {
      CP::THandle<CP::TG4HitContainer> g4Hits = (*h)->Get<CP::TG4HitContainer>(".");
      if (!g4Hits) {
        printf("truth/g4Hits object that is not a TG4HitContainer: %s ",(*h)->GetName());
        continue;
      }
      for (CP::TG4HitContainer::const_iterator h = g4Hits->begin();  h != g4Hits->end(); ++h) {
        const CP::TG4HitSegment* seg  = dynamic_cast<const CP::TG4HitSegment*>((*h));
        if (!seg) {
          printf("Not showing TG4Hit not castable as a TG4HitSegment!!! \n");
          continue;
        }

        CP::TCapHitSegment tseg;
        for( int icont = 0; icont< seg-> GetContributorCount(); ++icont) 
          tseg.Contributors.push_back(seg->GetContributor(icont));
        tseg.SetName((*h)->GetName());
        tseg.PrimaryId = seg->GetPrimaryTrajectoryId();
        tseg.EnergyDeposit = seg->GetEnergyDeposit();
        tseg.SecondaryDeposit = seg->GetSecondaryDeposit();
        tseg.TrackLength = seg->GetTrackLength();
        tseg.dedx = tseg.EnergyDeposit;
        if (tseg.TrackLength>0.01*unit::mm) tseg.dedx /= tseg.TrackLength;
        tseg.Start.SetXYZ(seg->GetStartX(),seg->GetStartY(),seg->GetStartZ());
        tseg.Stop.SetXYZ(seg->GetStopX(),seg->GetStopY(),seg->GetStopZ());
        tseg.StartT = seg->GetStartT();
        tseg.StopT = seg->GetStopT();

        
        int trajIndex = seg->GetContributor(0);
        if(fTrajectories) {
          CP::THandle<CP::TG4Trajectory> theTraj = fTrajectories->GetTrajectory(trajIndex);
          tseg.trajIndex = theTraj->GetTrackId();		
          /* these seem to be always the same but I am keeping the TrackId to be sure */
          //printf("  xxxxx %i =? %i xxxxx ",trajIndex,tseg.trajIndex);
          tseg.pdg = theTraj->GetPDGEncoding();
          tseg.SetTitle( theTraj->GetParticleName().c_str());
        }
        //tseg.print();
        fsegments.push_back(tseg);
      }
    }
    //printf(" size of fsegments %i \n",(int) fsegments.size());
  } // got truth hits 
  //printf(" here are the truth hits \n");
  //for(unsigned  ihit=0; ihit < fsegments.size(); ++ihit) fsegments[ihit].print(event.GetEventId());
  //now trajectories TrajectoriesModule
  // now we don't need the 
  //std::vector<CP::TG4TrajectoriesModule::TG4Trajectory> modTrajectories = fTG4TrajectoriesModule->fTrajectories;

  // trajectory number.  Iterating through the trajectory returns pairs with
  /// the first element being the trajectory number, and the second element
  /// being the TG4Trajectory object.
  
  CP::TG4TrajectoryContainer::iterator traj_it =  fTrajectories->begin();
  CP::TG4TrajectoryContainer::iterator traj_end = fTrajectories->end();
  CP::TG4Trajectory* traj;

  int ntraj=0;
  //printf(" filling trajectories \n ");
  int ntrajSegments=0;
  for (; traj_it != traj_end; ++traj_it) {
    ++ntraj;
    traj = &traj_it->second;
    //for(unsigned it =0; it<fTrajectories.size(); ++it) {
    CP::TCapTrajectory atraject;
    atraject.TrajId=traj->GetTrackId();		
    atraject.ParentId=traj->GetParentId();	
    atraject.PrimaryId=fTrajectories->GetPrimaryId(atraject.TrajId);		
    atraject.PDG=traj->GetPDGEncoding();
    atraject.Mass = CP::utils::FindPDGMass(traj->GetPDGEncoding());
    atraject.Charge = CP::utils::FindPDGCharge(traj->GetPDGEncoding());
    int ndrift=0; int ncryo=0; int nother=0;
    fillTruthPoints(&atraject,traj,ndrift,ncryo,nother);
    tcaptruth->NDrift=ndrift;
    tcaptruth->NCryo=ncryo;
    tcaptruth->NOther=nother;
    //add hits to this track 
    //atraject.print(ntraj);
    for(unsigned  iseg=0; iseg < fsegments.size(); ++iseg) { 
      //fsegments[iseg].print(iseg);
      if(  fsegments[iseg].trajIndex ==  atraject.TrajId ) atraject.segments.push_back(fsegments[iseg]);
    }
    //printf(" traj %i has %i segments \n ",ntraj,(int) atraject.segments.size());
    ntrajSegments += atraject.segments.size();
    //atraject.print(ntraj);
    //atraject.print(event.GetEventId());
    tcaptruth->trajectory.push_back(atraject);  
    //printf(" \n\n  blah blah size of pvertex %i trajectory % i fhits %i \n",(int) tcaptruth->pvertex.size() , (int) tcaptruth->trajectory.size(), 
      //  (int) fhits.size() );
  }
  printf(" \n  blah blah size of pverticies %i trajectories %i trajectory segments %i particles %i drift %i cryo %i other %i   \n",
      (int) tcaptruth->pvertex.size(),(int) tcaptruth->trajectory.size(),ntrajSegments , (int) tcaptruth->pparticle.size(),
      tcaptruth->NDrift,tcaptruth->NCryo ,tcaptruth->NOther );
  //tcaptruth->print();

}


void CP::captAna::Initialize(void){ 
  isFirstEvent=true;
  total_events=0;
  nWaveHist=0;
  for( int itype=0; itype<=UNKNOWN; ++itype) typeList[itype]=0;

  // list of class that inherit from TReconBase
  reconClassName.push_back(std::string("CP::TReconTrack"));
  reconClassName.push_back(std::string("CP::TReconShower"));
  reconClassName.push_back(std::string("CP::TReconCluster"));
  reconClassName.push_back(std::string("CP::TReconVertex"));
  reconClassName.push_back(std::string("CP::TReconPID"));
  printf(" initialized TReconBase classes :");
  for(unsigned iclass=0; iclass<reconClassName.size(); ++iclass) printf("  %s ",reconClassName[iclass].c_str());
  printf(" \n");

  gDirectory->cd("/");
  ntPmt = new TNtuple("ntAnaPmt"," pmt info ","ev:nhit:qsum");
  for(int jhist=0; jhist<MAXHIST; ++jhist) {
    hWaveSamples[jhist]=NULL;
    hDWaveSamples[jhist]=NULL;
    hElecResponse[jhist]=NULL;
    hFFTSamples[jhist]=NULL;
  }
  ilate=0;// non zero for charge injection data
  // new class based containers
  anaTree = new TTree("anaTree","anaTree");
  tcapevent = new CP::TCapEvent();
  tcaptruth = NULL;
  anaTree->Branch("event",&tcapevent);

  // this is old containers
  tdrift = new TTree("tdrift","tdrift");
  tdrift->Branch("rawWave",&rawWave,
      "iev/I:crate/I:fem/I:asic/I:nchan/I:plane/I:wire/I:baseMedian/F:baseSigma/F:sampleSigma/F:fitMean/F:fitMeanErr/F:fitSlope/F:fitSlopeErr/F:fitChi/F:npower/F");
  //maxadc=pow(2,12);//4096
  maxadc=9596;
  int maxSample=maxadc-1;
  //int halfadc=pow(2,11);
  hWaveFit = new TH1D("hWaveFit","wave to fit",maxSample,0,maxSample);
  hWaveFit->Sumw2();
  hClusterQ = new TH1D("ClusterQ"," cluster charge ",1000,0,1000000);
  int nSize = int(NSIZE);

  // FFT by FEM
  for(int ifem = 0; ifem<MAXFEM; ++ifem) {
    for(int jfem = 0; jfem<MAXASIC; ++jfem) {
      hFFT_FEM[ifem][jfem]=NULL;
      hFFT_DEC[ifem][jfem]=NULL;
    }
  }

  // initialze FFT
  fFFT = TVirtualFFT::FFT(1, &nSize, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nSize, "C2R M K");
  fFFTResponse= TVirtualFFT::FFT(1, &nSize, "R2C M K");
  fElectronicsResponse = new CP::TElectronicsResponse(nSize);      
  printf(" captAna with MAXHIST = %i \n",MAXHIST);
  gDirectory->cd("/");
  gDirectory->ls();
}

/******************************************/ 
// **** loop over events ***
/******************************************/ 
bool CP::captAna::operator () (CP::TEvent& event) {
  // clear the branches
  tcapevent->clear();
  if(tcaptruth) tcaptruth->clear();

  CP::TChannelInfo::Get().SetContext(event.GetContext());
  CP::TManager::Get().Geometry();

  CP::TChannelCalib channelCalib;
  int nSize = int(NSIZE);
  int ievent = event.GetEventId();
  tcapevent->EventId = event.GetEventId();
  tcapevent->RunId = event.GetRunId();
  ++total_events;


  // PMT geometry 
  if(isFirstEvent) {
    printf("GEOGEOGEOGEOGEOGEO  \n\n\n ");

   // Iterate through all of the defined geometry identifiers.
    for (CP::TGeomIdManager::GeomIdMap::const_iterator g
             = CP::TManager::Get().GeomId().GetGeomIdMap().begin();
         g != CP::TManager::Get().GeomId().GetGeomIdMap().end();
         ++g) {
        CP::TGeometryId geomId(g->first);
        if (!geomId.IsValid()) {
            std::cout << "Invalid geometry id: " << g->first << std::endl;
            continue;
        }
        if(geomId.GetName().find("Phot")!=string::npos) std::cout << (unsigned) g->first << " " << geomId.GetName() << std::endl;
    }
  }
    

  if(ievent%1==0) {
    printf(" ################ run %i event %i events read %u ################## \n",tcapevent->RunId,tcapevent->EventId,total_events);
    CaptLog("Event " << event.GetContext());
    CaptLog("Event " << event.GetContext());
    std::cout << " runid " << event.GetRunId()<< " ";
    std::cout << " timestamp " << event.GetContext().GetTimeStamp();
    std::cout << " nanosec " << event.GetContext().GetNanoseconds()<<std::endl;
    printf(" #################################################### \n");
  }
  tcapevent->timestamp=event.GetContext().GetTimeStamp();
  tcapevent->nanoseconds=event.GetContext().GetNanoseconds();
  //tcapevent->print();
  //std::cout << event.GetContext().GetSpill()<<std::endl
  if(ievent<2){
    printf(" data structure for this file : \n ");
    for (CP::TEvent::iterator v= event.begin(); v != event.end(); ++v) {
      printf(" TEvent iter: %s %s  \n",(*v)->GetName(),(*v)->GetTitle());
      // this gives more details
      (*v)->ls();
    }
    printf(" ----------------------  \n\n\n ");
  }


  // **** check for digits data vector ***
  if(isFirstEvent) printf(" checking for digits event %i \n",ievent);
  CP::THandle<CP::TDataVector> digits = event.Get<CP::TDataVector>("digits");
  if(!digits&&isFirstEvent) printf(" digtits is null event %i \n",ievent);
  if(isFirstEvent) {
    printf(" \n\n ---------- structure of digits data vector  : \n ");
    printf(" >>> size of digits %i  \n",int(digits->size()));
    for (CP::TDataVector::iterator dv= digits->begin(); dv != digits->end(); ++dv) {
      std::cout << "      ls from this iterator ";  (*dv)->ls();
      printf("  _____ TDataVector iter: %s %s  \n",(*dv)->GetName(),(*dv)->GetTitle());
    }
    printf(" ----------------------  \n\n\n ");
  }


  // **** check for digits data vector ***
  if(isFirstEvent) printf(" checking for fits event %i \n",ievent);
  CP::THandle<CP::TDataVector> dfits = event.Get<CP::TDataVector>("fits");
  if(!dfits&&isFirstEvent) printf(" digtits is null event %i \n",ievent);
  if(isFirstEvent) {
    printf(" \n\n ---------- structure of fits data vector  : \n ");
    printf(" >>> size of fits %i  \n",int(dfits->size()));
    for (CP::TDataVector::iterator df= dfits->begin(); df != dfits->end(); ++df) {
      std::cout << "      ls from this iterator ";  (*df)->ls();
      printf("  _____ TDataVector iter: %s %s  \n",(*df)->GetName(),(*df)->GetTitle());
    }
    printf(" ----------------------  \n\n\n ");
  }
  

  // look for pmt hits
/*
--------- structure of digits data vector  : 
  >>> size of digits 2  
      ls from this iterator CP::TDigitContainer(0x1af00d0):: pmt
 Signature: 4103335659
  _____ TDataVector iter: pmt Digit Pointers  
      ls from this iterator CP::TDigitContainer(0x1b41120):: drift-deconv
 Signature: 479930511
  _____ TDataVector iter: drift-deconv Digit Pointers  
 ----------------------  
*/
  CP::THandle<CP::TDigitContainer> pmtData = event.Get<CP::TDigitContainer>("digits/pmt");
  if(pmtData) {
    //printf(" size of pmt data %i  \n", (int) pmtData->size());
    for (CP::TDigitContainer::iterator tpmt = pmtData->begin(); tpmt != pmtData->end(); ++tpmt) {
      const CP::TPulseDigit *pdigi = dynamic_cast<const CP::TPulseDigit*>(*tpmt);
      if(!pdigi)  { printf(" no pdigi found !!!!! "); continue; }
      //pdigi->ls();
    }
  }

  /** this is the raw data **/
  CP::THandle<CP::TDigitContainer> drift = event.Get<CP::TDigitContainer>("digits/drift");
  if(!drift&&isFirstEvent) printf(" drift is null event %i \n",ievent);
  if(drift) {
    if(isFirstEvent) printf(" size of drift selection %i  \n",(int) drift->size());
    //for (std::size_t d = 0; d < drift->size(); ++d) {}
    int num_channels=0;
    for (CP::TDigitContainer::iterator td = drift->begin(); td != drift->end(); ++td) {
      //const CP::TPulseDigit* digit= dynamic_cast<const CP::TPulseDigit*>((*drift)[d]);
      const CP::TPulseDigit *digit = dynamic_cast<const CP::TPulseDigit*>(*td);
      if(!digit) {
        printf(" no digit found !!!!! ");
        continue;
      }
      //digit->ls();
      //printf(" sample count %i \n",(int)digit->GetSampleCount());
      // loop over sample counts and histogram
      CP::TChannelId chid = digit->GetChannelId();
      double timeStep = channelCalib.GetTimeConstant(chid);
      // unit::nanosecond is one
      if(isFirstEvent&&num_channels==0) cout << chid.AsString() << " timestep is  "
        << timeStep*unit::nanosecond << " nsec  " << timeStep/unit::microsecond << " micro-sec " << endl; 
      ++num_channels;
      CP::TTPCChannelId tpcid; 
      if(chid.IsTPCChannel())  tpcid = CP::TTPCChannelId(chid);
      else continue;
      vector<unsigned short> samples = digit->GetSamples();

      // calculate baseline 
      //< clusterCalib.deconvolution.baselineCut = 5 >

      // make a vector to sort the samples
      std::vector<double> diff;
      diff.resize(digit->GetSampleCount());
      // Find the sample median and it's "sigma".
      for (std::size_t i=1; i<digit->GetSampleCount(); ++i) {
        hWaveFit->SetBinContent(i,double(digit->GetSample(i)));  // for fitting
        hWaveFit->SetBinError(i,1.0); // one adc count is the error
        diff[i] = digit->GetSample(i);
      }
      std::sort(diff.begin(), diff.end());

      double baselineMedian = diff[0.5*digit->GetSampleCount()];
      double baselineSigma = diff[0.16*digit->GetSampleCount()];
      baselineSigma = std::abs(baselineSigma-baselineMedian);

      //fit the wave to a line

      TF1 *fline = new TF1("fline", "pol1"); // a line 
      fline->SetParameter(0,baselineMedian);
      fline->SetParameter(1,0);
      TFitResultPtr fitptr = hWaveFit->Fit("fline","QS"); // M for minuit L for likelihood Q for quiet N is to not store result E minos errors
      double fitMean    = fitptr->Parameter(0);
      double fitMeanErr = fitptr->ParError(0);
      double fitSlope   = fitptr->Parameter(1);
      double fitSlopeErr= fitptr->ParError(1);
      double fitChi = fitptr->Chi2()/fitptr->Ndf();

      // The minimum baseline fluctuation is 1 electron charge.
      if (baselineSigma < 1) baselineSigma = 1.0;

      // Define a maximum separation between the sample and the median sample.
      // If it's more than this, then the sample is not baseline.  This is one
      // sided.  This is looking at the overall fluctation of the baseline and
      // prevents large plateaus from being cut.
      //double baselineCut = baselineMedian + fBaselineCut*baselineSigma;

      // Find the median sample to sample difference.  Regions where the samples
      // stay withing a small difference don't have a "feature of interest".
      for (std::size_t i=1; i<digit->GetSampleCount(); ++i) {
        double delta = std::abs(digit->GetSample(i) - digit->GetSample(i-1));
        diff[i] = delta;
      }
      std::sort(diff.begin(), diff.end());

      // The 52 percentile is a of the sample-to-sample differences is a good
      // estimate of the distribution sigma.  The 52% "magic" number comes
      // looking at the differences of two samples (that gives a sqrt(2)) and
      // wanting the RMS (which is the 68%).  This gives 48% (which is
      // 68%/sqrt(2)).  Then because of the ordering after the sort, the bin we
      // look at is 1.0-52%.  The minimum allowable fluctuation is 1 electron
      // charge.
      double sampleSigma = diff[0.52*digit->GetSampleCount()];
      if (sampleSigma < 1.0) sampleSigma = 1.0;

      //
      // make list of histograms to save
      if(nWaveHist<MAXHIST) {
        if(tpcid.GetCrate()==2&&tpcid.GetFEM()==11&&ievent==9) {
          tpcid_list[nWaveHist]=tpcid; 
          event_list[nWaveHist]=ievent;
          ++nWaveHist;
        }
      }

      // fill waveforms in list 
      //char hname[80];
      char htitle[80];
      int maxSample=maxadc-1;


      CP::TGeometryId geoid = CP::TChannelInfo::Get().GetGeometry(chid);
      int wireNum = int(CP::GeomId::Captain::GetWireNumber(geoid) + 0.5);
      

      for(int jhist=0; jhist<nWaveHist; ++jhist) {
         if(!hWaveSamples[jhist]) {
           std::string stpcid=ChIdString(tpcid_list[jhist]);
           sprintf(htitle,"rWave_%s_Wire_%i_Ev_%i",stpcid.c_str(),wireNum,event_list[jhist]);
           hWaveSamples[jhist] = new TH1D(htitle,htitle,maxSample,0,maxSample);
           sprintf(htitle,"dWave_%s_Wire_%i_Ev_%i",stpcid.c_str(),wireNum,event_list[jhist]);
           hDWaveSamples[jhist] = new TH1D(htitle,htitle,maxSample,0,maxSample);
           printf(" created histograms %s %s \n",hWaveSamples[jhist]->GetName(),hDWaveSamples[jhist]->GetName());
         }
         if(hWaveSamples[jhist]&&tpcid_list[jhist]==tpcid&&event_list[jhist]==ievent) 
           for(unsigned is =0; is<samples.size(); ++is) hWaveSamples[jhist]->SetBinContent(is+1,double(samples[is]));
      }

      // **** loop over samples , TDC times for FFT 
      //  fft for this digit 
      //printf(" setting fft points \n\n");
      for(unsigned is =0; is<samples.size(); ++is) fFFT->SetPoint(is,samples[is]); // FFT of signal
      fFFT->Transform();


      CP::TChannelCalib channelCalib;
      // look at Electronic Response 
      vector<double> fResponse(NSIZE);
      // Fill the response function.  This explicitly normalizes so that the
      // rawWave shaping is amplitude conserving (the rawWave shaping for a
      // "delta-function" sample doesn't change the rawWave height).
      double normalization = 0.0;
      for (unsigned i=0; i<fResponse.size(); ++i) {
        //double arg = channelCalib.GetTimeConstant(chid)*(1.0*i+0.5);
        //double v = channelCalib.GetPulseShape(chid,arg);
        //fResponse[i] = v;
        fResponse[i] = fElectronicsResponse->GetResponse(i);
        normalization = std::max(normalization,fResponse[i]);
      }
      assert(normalization > 1E-20);
      // normalize 
      for (std::vector<double>::iterator r = fResponse.begin(); r != fResponse.end(); ++r) {
        (*r) /= normalization;
      }

      // need fft of response function for noise calculation 
      for(unsigned is =0; is<fResponse.size(); ++is) fFFTResponse->SetPoint(is,fResponse[is]);      // **** loop over samples , TDC times for FFT 
      fFFTResponse->Transform();

      // fill samples FFT histogram && elec response in time domain
      for(int jhist=0; jhist<MAXHIST; ++jhist) {
        if(tpcid == tpcid_list[jhist] ) {
          std::string chanName = ChIdString(tpcid);
          if(!hFFTSamples[jhist]){ 
            hFFTSamples[jhist] = new TH1D((chanName+std::string("_fft")).c_str(),
                ("FFT for " + chid.AsString()).c_str(),
                nSize/2,0,(0.5/timeStep)/unit::hertz);
            printf(" created %i %s %s \n",jhist,hFFTSamples[jhist]->GetName(),hFFTSamples[jhist]->GetTitle());
          }
          // skip first bin which is pedestal
          for (int i = 1; i<nSize/2; ++i) {
            double rl, im;
            fFFT->GetPointComplex(i,rl,im);
            std::complex<double> c(rl,im);
            hFFTSamples[jhist]->SetBinContent(i+1,hFFTSamples[jhist]->GetBinContent(i+1)+std::abs(c));
          } 

          // elec response in the time domain
          if(!hElecResponse[jhist]){
            hElecResponse[jhist] = new TH1D((chanName+std::string("_ElecResponse")).c_str(),
                ("elec response " + chid.AsString()).c_str(),nSize,0,nSize);
            for (unsigned i = 1; i<fResponse.size(); ++i) hElecResponse[jhist]->SetBinContent(i+1,fResponse[i]);
            printf(" created %i %s %s \n",jhist,hElecResponse[jhist]->GetName(),hElecResponse[jhist]->GetTitle());
          }
        }
      }

      // ***************** moise filter statistics *********************
      vector<double> fWork(NSIZE); // abs value of signal
      for (std::size_t i = 0; i<NSIZE; ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        fWork[i] = std::abs(c);
      }

      vector<double> fRWork(NSIZE);  // abs value of response
      for (std::size_t i = 0; i<NSIZE; ++i) {
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        fRWork[i] = std::abs(c);
      }

      double maxResponse = 0.0;
      double minResponse = 1E+6;
      for (std::size_t i = 0; i<NSIZE; ++i) {
        maxResponse = std::max(maxResponse,fRWork[i]);
        minResponse = std::min(minResponse,fRWork[i]);
      }

      // Find the Gaussian noise estimate.
      double maxSignal = 0.0;
      double maxSignalWeight = 0.0;
      double minSignal = 0.0;
      double minSignalWeight = 0.0;
      for (unsigned i=0; i<fWork.size(); ++i) {
        double w = fRWork[i]/maxResponse;
        w = w*w;
        double p = fWork[i]*fWork[i]/fWork.size()/fWork.size();
        maxSignalWeight += w;
        maxSignal += w*p;
        w = 1.0-fRWork[i]/maxResponse;
        w *= w;
        minSignalWeight += w;
        minSignal += w*p;
      }
      maxSignal = maxSignal/maxSignalWeight;
      minSignal = minSignal/minSignalWeight;
      maxSignal = std::sqrt(maxSignal*fWork.size()/2.0);
      minSignal = std::sqrt(minSignal*fWork.size()/2.0);

      // Estimate the noise power relative to the averaged signal.
      double npower = 0.0;
      if (minSignal < maxSignal) {
        npower = maxSignal*minResponse-maxResponse*minSignal;
        npower /= (minSignal-maxSignal);
      }

      // fill by FEM/ASIC
      int nChannel = (int) tpcid.GetChannel();
      int nFEM  =tpcid.GetFEM();
      int nASIC=nChannel>>4;

      // get geo plane and wire  number.
      //CP::TGeometryId geoid = CP::TChannelInfo::Get().GetGeometry(chid);
      int plane = CP::GeomId::Captain::GetWirePlane(geoid);
      double wire = CP::GeomId::Captain::GetWireNumber(geoid) + 0.5;

      char hname[80];
      for (int i = 1; i<nSize/2; ++i) {
        if(!hFFT_FEM[nFEM][nASIC]) {
          sprintf(hname,"FFT_FEM%02iASIC%02i",nFEM,nASIC);
          sprintf(htitle,"FFT for FEM %02i ASIC %02i",nFEM,nASIC);
          //printf(" %s %s \n",hname,htitle);
          hFFT_FEM[nFEM][nASIC] = new TH1D(hname,htitle,nSize/2,0,(0.5/timeStep)/unit::hertz);
        }
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        hFFT_FEM[nFEM][nASIC]->SetBinContent(i+1, hFFT_FEM[nFEM][nASIC]->GetBinContent(i+1)+std::abs(c));
      } 

      rawWave.iev=ievent;
      rawWave.crate=tpcid.GetCrate();
      rawWave.wire = int(wire);
      rawWave.plane = plane;
      rawWave.fem= nFEM;
      rawWave.asic=nASIC;
      rawWave.nchan= tpcid.GetChannel();
      rawWave.baseMedian = baselineMedian;
      rawWave.baseSigma  = baselineSigma;
      rawWave.sampleSigma  = sampleSigma;
      rawWave.fitMean= fitMean;
      rawWave.fitMeanErr= fitMeanErr;
      rawWave.fitSlope= fitSlope;
      rawWave.fitSlopeErr = fitSlopeErr;
      rawWave.fitChi = fitChi;
      rawWave.npower = float(npower);
      tdrift->Fill();
    }  // raw digit container loop 
  }// if drift 

  //for(int jhist=0; jhist<MAXHIST; ++jhist) printf(" %s \n",tpcid_list[jhist].AsString().c_str());
  CP::THandle<CP::TDigitContainer> dcdrift = event.Get<CP::TDigitContainer>("digits/drift-deconv");
  if(!dcdrift&&isFirstEvent) printf(" dcdrift is null event %i \n",ievent);
  if(dcdrift) {
    if(isFirstEvent) printf(" size of dcdrift selection %i  \n",(int) dcdrift->size());
    for (CP::TDigitContainer::iterator tdc = dcdrift->begin(); tdc != dcdrift->end(); ++tdc) {
      CP::TCalibPulseDigit *dcdigit = (CP::TCalibPulseDigit *) (*tdc);
      if(!dcdigit) continue;
      //dcdigit->ls();
      CP::TChannelId chID = dcdigit->GetChannelId();
      CP::TTPCChannelId tpcID; 
      if(chID.IsTPCChannel())  tpcID = CP::TTPCChannelId(chID);
      //else continue;// if not a TPCChannel 
      //int icrate = tpcID.GetCrate();

      vector<float> samples = dcdigit->GetSamples();
      // loop over all samples for FFT
      for(unsigned is =0; is<samples.size(); ++is)  fFFT->SetPoint(is,samples[is]);
      fFFT->Transform();
      // fill by FEM/ASIC
      int nChannel =tpcID.GetChannel();
      int nFEM  =tpcID.GetFEM();
      int nASIC=nChannel>>4; 
      char hname[80];
      char htitle[80];
      double timeStep = channelCalib.GetTimeConstant(chID);

      for (int i = 1; i<nSize/2; ++i) {
        if(!hFFT_DEC[nFEM][nASIC]) {
          sprintf(hname,"FFT_DEC_FEM%02iASIC%02i",nFEM,nASIC);
          sprintf(htitle,"deconvolved FFT for FEM %02i ASIC %02i",nFEM,nASIC);
          //printf(" %s %s \n",hname,htitle);
          hFFT_DEC[nFEM][nASIC] = new TH1D(hname,htitle,nSize/2,0,(0.5/timeStep)/unit::hertz);
        }
        double rl, im;
        fFFT->GetPointComplex(i,rl,im);
        std::complex<double> c(rl,im);
        hFFT_DEC[nFEM][nASIC]->SetBinContent(i+1, hFFT_DEC[nFEM][nASIC]->GetBinContent(i+1)+std::abs(c));
      } 


      // fill waveforms in list 
      for(int jhist=0; jhist<MAXHIST; ++jhist) {
         if(hDWaveSamples[jhist]&&tpcid_list[jhist]==tpcID&&event_list[jhist]==ievent) 
          for(unsigned is =0; is<samples.size(); ++is) hDWaveSamples[jhist]->SetBinContent(is+1,double(samples[is]));
      }
    } // dcdigit
  }// if dcdrift

  CP::THandle<CP::THitSelection>  hits = event.GetHits("drift");  // gets hits/drift
  //CP::THandle<CP::THitSelection> hits = event.GetHits("hits");
  if(!hits&&isFirstEvent) printf(" hits is null event %i \n",ievent);
  if(hits) {
    if(isFirstEvent) printf(" size of hits selection %i  \n",(int) hits->size());
    TH1D* hhitwave = NULL;
    TH1D* hrawwave = NULL;
    TH1D* hdecwave = NULL;

    std::string current_stpcid;
    int hit_count=0;
    CP::TChannelId currentID(0);
    double channel_charge=0;
    char hname[120];
    char htitle[120];

    int num_hit=0;
    for (CP::THitSelection::iterator hhit = hits->begin(); hhit != hits->end(); ++hhit) {
      CP::THandle<CP::TPulseHit> theHit = *hhit;


      CP::TChannelId chid(theHit->GetChannelId());
      if(!chid.IsValid()) {
        printf("NOT a valid fucking TChannelID %s \n",chid.AsString().c_str());
        continue;
      }

      // get geo plane and wire  number.
      CP::TGeometryId geoid = CP::TChannelInfo::Get().GetGeometry(chid);
      //int plane = CP::GeomId::Captain::GetWirePlane(geoid);
      //double wire = CP::GeomId::Captain::GetWireNumber(geoid) + 0.5;


      CP::TTPCChannelId tpcid(chid);
      std::string stpcid=ChIdString(chid);
      // fill channel info for the capHit
      //int nChannel =tpcid.GetChannel();
      //int nASIC= nChannel>>4; 

      // check if new channel
      if(!(chid==currentID)) {
        sprintf(htitle,"hitWave_%s_Ev_%i charge %E hits %i",current_stpcid.c_str(),ievent,channel_charge,hit_count);
        if(hhitwave) hhitwave->SetTitle(htitle);
        currentID=chid;
        current_stpcid=stpcid;
        hhitwave=NULL;
        hrawwave=NULL;
        hdecwave=NULL;
        channel_charge=0;
        hit_count=0;
      }
      //double timeStep = channelCalib.GetTimeConstant(chid);

      // take nominal pedestal for now.
      double pedestal = channelCalib.GetGainConstant(chid,0);
      double gain = channelCalib.GetGainConstant(chid,1);
      double slope = channelCalib.GetDigitizerConstant(chid,1);

      //printf("raw wave ped %E gain %E slope %E \n",pedestal,gain,slope);


      CP::TDigit *rawDigit=NULL; 
      if(drift) rawDigit = GetDigit( drift, chid );
      else printf(" aint got no drift \n");


      CP::TDigit *decDigit=NULL; 
      if(dcdrift) decDigit = GetDigit( dcdrift, chid );
      //else printf(" aint got no deconvolved drift \n");


      // calculate the noise and baseline
      double fSampleSigma=0;
      if(decDigit) {
        CP::TCalibPulseDigit *dcdigit = (CP::TCalibPulseDigit *) (decDigit);
        vector<float> diff = dcdigit->GetSamples();
        std::sort(diff.begin(), diff.end());
        // The 52 percentile is a of the sample-to-sample differences is a good -- see TPulseDeconvolution...
        int sigma52 = int(0.52*dcdigit->GetSampleCount());
        fSampleSigma = diff[sigma52];
        if (fSampleSigma < 1.0) fSampleSigma = 1.0;
      }


      // get corresponding raw hit
      channel_charge += theHit->GetCharge();
      ++hit_count;
      if(ievent<20&&rawDigit&&decDigit) {
        CP::TPulseDigit *rawPulse =   dynamic_cast<CP::TPulseDigit *>(rawDigit);
        CP::TCalibPulseDigit *decPulse =   dynamic_cast<CP::TCalibPulseDigit *>(decDigit);
        //printf("  found this rawWave \n");
        //rawPulse->ls();

        vector<unsigned short> rawSamples = rawPulse->GetSamples();
        vector<float> decSamples = decPulse->GetSamples();
        int maxSample = rawPulse->GetSampleCount();

        //printf("  found in event %i this rawWave  %s  charge %E hit %i  \n",
        //ievent,tpcid.AsString().c_str(),theHit->GetCharge(),hit_count);
        if(!hhitwave) {
          sprintf(hname,"hitWave_%s_Ev_%i",stpcid.c_str(),ievent);
          sprintf(htitle,"hitWave_%s_Ev_%i charge %E hits %i",stpcid.c_str(),ievent,channel_charge,hit_count);
          hhitwave = new TH1D(hname,htitle,maxSample,0,maxSample);
        }
        sprintf(htitle,"rawWave_%s_Ev_%i",stpcid.c_str(),ievent);
        if(!hrawwave) hrawwave = new TH1D(htitle,htitle,maxSample,0,maxSample);
        sprintf(htitle,"decWave_%s_Ev_%i",stpcid.c_str(),ievent);
        if(!hdecwave) hdecwave = new TH1D(htitle,htitle,maxSample,0,maxSample);

        //double hitTime = theHit->GetTime();
        //double hitRMS = theHit->GetTimeRMS();

        // Find the time per sample in the digit.
        double digitStep = decPulse->GetLastSample()-decPulse->GetFirstSample();
        digitStep  /= decPulse->GetSampleCount();

        //printf(" digitStep %f timeStep %f first sample %f \n",digitStep,timeStep,decPulse->GetFirstSample()); 
        // //these are the same!! = 500.0  first sample is -1600000.000000 

        for(unsigned is =0; is<rawSamples.size(); ++is) {
          double stime = digitStep*(double(is)+0.5)+decPulse->GetFirstSample();
          if( stime > theHit->GetTimeStart() && stime < theHit->GetTimeStop() ) {
            hhitwave->SetBinContent(is+1,hhitwave->GetBinContent(is+1)+double(decSamples[is]));
          }
        }

        //printf(" <<<<<< %s %i event %i  \n",hWaveSamples[jhist]->GetTitle(),jhist,event_list[jhist]);
        for(unsigned is =0; is<rawSamples.size(); ++is){
          hrawwave->SetBinContent(is+1,(double(rawSamples[is])-pedestal)/gain/slope);
          hdecwave->SetBinContent(is+1,double(decSamples[is]));
        }
      } //else 
      //printf(" digit not found %s \n",chid.AsString().c_str());

      hClusterQ->Fill( (*hhit)->GetCharge());
      CP::TCapHit capHit;
      fillCapHit(&capHit,theHit,fSampleSigma);
      //capHit.print(num_hit);
      tcapevent->hits.push_back(capHit);
      ++num_hit;
    }
  }

  // look for PMT hits
  CP::THandle<CP::THitSelection>  phits = event.GetHits("pmt");  // gets hits/drift
  if(!phits) {
    printf(" PMT hits is null event %i \n",ievent);
  } else  { 
    //printf(" \n \n blahblahblah size of PMT hits selection %i  \n",(int) phits->size());
    //phits->ls();

    int npmtHits =0;
    float qsumHits=0;
    for (CP::THitSelection::iterator piter = phits->begin(); piter != phits->end(); ++piter) {
      CP::TCapPMTHit capPHit;
      capPHit.geom = (*piter)->GetGeomId(0).GetName();
      capPHit.time = (*piter)->GetTime();
      capPHit.timeUnc = (*piter)->GetTimeUncertainty();
      capPHit.timeRms = (*piter)->GetTimeRMS();
      capPHit.charge = (*piter)->GetCharge();
      ++npmtHits;
      qsumHits += (*piter)->GetCharge();
      capPHit.chargeUnc =(*piter)->GetChargeUncertainty();
      capPHit.position =(*piter)->GetPosition();
      capPHit. positionUnc =(*piter)->GetUncertainty();
      capPHit.positionRms =(*piter)->GetRMS();
      capPHit.tstart =(*piter)->GetTimeStart();
      capPHit.tstop =(*piter)->GetTimeStop();
      capPHit.nsamples =(*piter)->GetTimeSamples();
      for(int isample =0; isample< capPHit.nsamples ;  ++isample ) 
        capPHit.qsamples.push_back( (*piter)->GetTimeSample(isample) );		

      //capPHit.print( tcapevent->phits.size() ); 

      tcapevent->phits.push_back(capPHit);

      // check
    }
    ntPmt->Fill(float(event.GetEventId()),float(npmtHits),qsumHits);
  }// if(phits)

  // look for fits
  CP::THandle<CP::TDataVector> fitptr =event.Get<CP::TDataVector>("~/fits/TCaptainRecon");
  if(fitptr&&isFirstEvent) printf(" xxxxxx  fits is FOUND full name is  %s \n", fitptr->GetFullName().Data());
  if(!fitptr&&isFirstEvent) printf(" xxxxxx TCaptainRecon not found!  \n");

  if(fitptr) {
    if(isFirstEvent) printf("  size of fits %i  \n",(int) fitptr->size());
    /*
       name  TCluster3D title An Algorithm Result  isAlgo 1 
       name  TClusterSlice title An Algorithm Result  isAlgo 1 
       name  TMinimalSpanningTrack title An Algorithm Result  isAlgo 1 
       name  TSplitTracks title An Algorithm Result  isAlgo 1 
       name  TMergeTracks title An Algorithm Result  isAlgo 1 
       name  TDisassociateHits title An Algorithm Result  isAlgo 1 
       name  TCombineOverlaps title An Algorithm Result  isAlgo 1 
       name  hits title Link To  isAlgo 0 
       name  unused title Hit Handles  isAlgo 0 
       name  used title Hit Handles  isAlgo 0 
       name  results title Link To  isAlgo 0 
       name  final title Recon Object Container  isAlgo 0 
       */
    std::vector<string> algoName;
    CP::TReconObjectContainer* fresult=NULL; //defined in TReconBase.hxx
    for (CP::TAlgorithmResult::iterator tv = fitptr->begin(); tv != fitptr->end(); ++tv) {
      TObject *obj = *tv;
      bool isAlgo = obj->InheritsFrom("CP::TAlgoritmResult");
      //printf(" TAlogrithmResult name  %s title %s  ? %i \n ",(*tv)->GetName(),(*tv)->GetTitle(),isAlgo);
      std::string aname((*tv)->GetName());
      if(strcmp(aname.c_str(),"final")==0) fresult = (CP::TReconObjectContainer*) (*tv);
      //(*tv)->ls();
      if(isAlgo) {
        algoName.push_back(aname);
        CP::TAlgorithmResult*  result = (CP::TAlgorithmResult*) (*tv);
        if(isFirstEvent) printf(" xxxx  name  %s title %s  \n",result->GetName(),result->GetTitle());
        //result->GetAlgorithmTag().ls();//
      }//
      //printf("\n\n");
    }
    // class CP::TReconObjectContainer : public TDatum, public std::vector< CP::THandle<CP::TReconBase> > {
    if(fresult) {
       if(total_events%10==0) printf("run %i event %i ................final results container size %i \n",event.GetRunId(),event.GetEventId(), (int) fresult->size());
      for (CP::TReconObjectContainer::iterator rv = fresult->begin(); rv != fresult->end(); ++rv) {
        reconType rtype=UNKNOWN;
        for(unsigned iclass=0; iclass<reconClassName.size(); ++iclass) { 
          if( (*rv)->InheritsFrom(reconClassName[iclass].c_str()) ) { 
            rtype=reconType(iclass);   
            ++typeList[rtype];      
             if(total_events%10==0) printf(" iclass %i type %i TReconObject %s \n",iclass, rtype,reconClassName[iclass].c_str());
             if(total_events%10==0) printTypeList();
          }
        }
        if(rtype==UNKNOWN) printf(" CAPTANA does not know of this TReconObject!! ");
        CP::THandle<CP::THitSelection> fhitselect = (*rv)->GetHits();
        if(rtype==RTRACK){
          //printf(" final result track   %s title %s hits %i  \n",(*rv)->GetName(),(*rv)->GetTitle(),(int) fhits->size() );
          CP::TCapTrack capTrack;
          CP::THandle<CP::TReconTrack> rtrack = *rv;
          //printf(" \n\n\n blah blah blah \n\n\n");
          //rtrack.ls();
          fillCapTrack(&capTrack,rtrack);
          //capTrack.print(1);
          tcapevent->tracks.push_back(capTrack);
        } else if(rtype==RCLUSTER) {
          CP::TCapCluster capCluster;
          CP::THandle<CP::TReconCluster> rcluster = *rv;
          fillCapCluster(&capCluster,rcluster);
          tcapevent->clusters.push_back(capCluster);
        } else if(rtype==RSHOWER) {
          CP::TCapShower capShower;
          CP::THandle<CP::TReconShower> rshower = *rv;
          fillCapShower(&capShower,rshower);
          tcapevent->showers.push_back(capShower);
        } else {
          printf(" NNNNNNNNNNNNNNNNNNNNNN  final result-- i dont have this type yet ?  %s title %s hits %i  \n",
              (*rv)->GetName(),(*rv)->GetTitle(),(int) fhitselect->size() );
        }
        //(*rv)->ls();
      }
      //fresult->ls();
      // these seem to be TReconCluster, TReconTrack
      // member of base class CP::THandle<CP::THitSelection> GetHits() 
      // get hits and save for this event


    }
    //printf(" alog names %i \n",(int) algoName.size());
    //for (CP::TReconObjectContainer::iterator tvfit = fits->begin(); fits != fits->end(); ++hhit) { }
    //CP::THandle<CP::THitSelection> rhits = fits->GetHits();
  }



  // look for truth 
  CP::THandle<CP::TDataVector> truthptr = event.Get<CP::TDataVector>("truth");
  if(truthptr) {
    if(isFirstEvent) printf(" XXXXX GOT THE TRUTH XXXXX size of truth %i  \n",(int) truthptr->size());
    if(isFirstEvent) {
      InitTruth();
      fTG4TrajectoriesModule->ProcessFirstEvent(event);
      fTG4VerticiesModule->ProcessFirstEvent(event);
    } 
    fTG4TrajectoriesModule->Process(event);
    fTG4VerticiesModule->Process(event);
    fillTruth(event);
  }
  // tree filled on every event
  if(truthptr&&tcaptruth->NDrift>0) {
    anaTree->Fill(); 
    if(total_events%10==0) printf(" \n\n  woof woof event %i  size of pverticies %i trajectories %i  particles %i drift %i cryo %i other %i anaTree size %i \n",
        total_events,(int) tcaptruth->pvertex.size(),(int) tcaptruth->trajectory.size(), (int) tcaptruth->pparticle.size(),
        tcaptruth->NDrift,tcaptruth->NCryo ,tcaptruth->NOther ,(int) anaTree->GetEntriesFast() );
  } else if(!truthptr){ 
    anaTree->Fill();
    if(total_events%10==0) printf(" \n\n  woof woof event %i  anaTree size %i \n",total_events, (int) anaTree->GetEntriesFast() );
  }

   
  isFirstEvent=false; 
  return false;
} // end of event routine


/******************************************/ 
void CP::captAna::Finalize(CP::TRootOutput * const output) {

  printf(" events processed %u\n",total_events);
  printTypeList();
  if(output) {
    TDirectory* truthDir = (TDirectory *) output->FindObject("TruthDir");
    if(truthDir) {
      printf(" Deleting redundant TruthDir \n");
      truthDir->Delete();
    }
    TTree* captainTree = (TTree *) output->FindObject("captainEventTree");
    if(captainTree) {
      printf(" Deleting captainEventTree \n");
      captainTree->Delete();
    }
    
    printf(" anaTree %i output events written %i \n", (int) anaTree->GetEntriesFast(), output->GetEventsWritten() ); 
  } else 
    printf(" output not written as no output file specified \n");
}


int main(int argc, char **argv) {
  cout << " executing " << argv[0] << endl;
  for(int i=1; i<argc; ++i) printf(" %i %s ",i,argv[i]);
  if(argc>2) printf(" output file %s \n",argv[3]);
  CP::captAna userCode;
  std::string inputFile(argv[argc-1]);
  cout << " geometry from file " << inputFile << endl;
  CP::TManager::Get().SetGeometryOverride(inputFile.c_str());  

  CP::eventLoop(argc,argv,userCode);
}

