/**
 *   \file
 *   implementation of class FlowQA
 *   \author S. Johnson
 *   \date 12/9/2012
 *   \author extended by Anar Rustamov
 */

#include "FlowQA.h"
#include <fwk/CentralConfig.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/UTCDateTime.h>
#include <utl/WithUnit.h>
#include <evt/Event.h>
#include <evt/RecEvent.h>
// #include <evt/SimEvent.h>
#include <evt/rec/RecEventConst.h>
#include <det/TriggerConst.h>
#include <evt/EventHeader.h>
#include <utl/TimeStamp.h>
#include <det/Detector.h>
#include <det/MagneticFieldTracker.h>
#include <det/TPC.h>
#include <det/PSD.h>
#include <evt/rec/BPDPlane.h>
#include <evt/rec/RecEventConst.h>
#include <det/BPDConst.h>

#include <TFile.h>
#include <TTree.h> //vblinov
#include <TMath.h> //vblinov
#include <TClonesArray.h> //vblinov
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TDirectory.h>

#include <det/TPCConst.h>//vblinov
//#include <my_point.h>

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace fwk;
using namespace evt;
using namespace utl;
using namespace std;
using namespace evt;
using namespace evt::rec;
using namespace evt::sim;

namespace FlowQA {
  
/// Init function for FlowQA
VModule::EResultFlag FlowQA::Init()
{
    CentralConfig& cc = CentralConfig::GetInstance();
    
    Branch topBranch = cc.GetTopBranch("FlowQA");
    InitVerbosity(topBranch);
    
    // Get the QA file name for output.
    topBranch.GetChild("qaFileName").GetData(fQAFileName);
    
    //Counters
    fNEvents = 0;
    fNT1 = 0;
    fNT2 = 0;
    fNT3 = 0;
    
    fNGoodBPD = 0;
    fNHasMainVertex = 0;
    fNMainVertexPosition = 0;
    fNHasFittedVertex = 0;
    fNVertexTracks = 0;
    
    fNTracks = 0;
    fNTrackStatus = 0;
    fNVTPCClusters = 0;
    fNTotClusters = 0;
    fNClusterRatio = 0;
    fNImpactParameter = 0;
    fNHasTOF = 0;
    fUseVertexTracks = 0;
    
    timeUsed = kFALSE;
    timeStart = 0;
    eBeam = 1.;
    plotDriftVel = 0;
    trigBit.reset();
    for(int i = 0 ; i < 5; ++i)
      dVel[i] = 0;
    
    
    // cut values - can be changed in the .xml
    
    topBranch.GetChild("eventCutsOn").GetData(fEventCutsOn);
    topBranch.GetChild("loopThroughMainVertex").GetData(fLoopThroughMainVertex);
    topBranch.GetChild("useVertexTracks").GetData(fUseVertexTracks);
    
    if (fEventCutsOn)
      topBranch.GetChild("nMinVertexTracks").GetData(fNMinVertexTracks);
    else
      fNMinVertexTracks = 0;
    
    topBranch.GetChild("nMinVTPCClusters").GetData(fNMinVTPCClusters);
    topBranch.GetChild("nMinTotalClusters").GetData(fNMinTotalClusters);
    topBranch.GetChild("minPointRatio").GetData(fMinPointRatio);
    topBranch.GetChild("maxImpactX").GetData(fMaxImpactX);
    topBranch.GetChild("maxImpactY").GetData(fMaxImpactY);
    topBranch.GetChild("zMin").GetData(fZMin);
    topBranch.GetChild("zMax").GetData(fZMax);
    topBranch.GetChild("eBeam").GetData(eBeam);
    topBranch.GetChild("plotDriftVel").GetData(plotDriftVel);
    
    TDirectory * tmpDir = (TDirectory*)gROOT->CurrentDirectory();
    outputFile =  new TFile(fQAFileName.c_str(), "RECREATE");    
    outputFile->cd();
    fTreeQA = new TTree("fTreeQA", "Output tree for QA");
//     fTreeQA->SetMaxTreeSize(90000000);
    Init_fTreeQA();
    tmpDir->cd();
    return eSuccess;
  }
  
  /// Process function for FlowQA
  VModule::EResultFlag
  FlowQA::Process(Event& event, const AttributeMap& /*attr*/)
  {
//       std::cout << "ClearEvent" << std::endl;
    ClearEvent();
//       std::cout << "Calculate_QA" << std::endl;
    Calculate_QA(event);
    return eSuccess;
  }
  
  
  /// Finish function for FlowQA
  VModule::EResultFlag
  FlowQA::Finish()
  {
    outputFile->cd();
    fTreeQA -> Write();
//     outputFile->Close();
    TDirectory * tmpDir = (TDirectory*)gROOT->CurrentDirectory();
    tmpDir -> cd();
    return eSuccess;
  }
  
  void FlowQA::Calculate_QA(Event& event)
  {
      cout << "FlowQA::Calculate_QA" << endl;
      //Event
      const EventHeader& evtHeader = event.GetEventHeader();
      Event_Timestamp = evtHeader.GetTime().GetGPSSecond();
      Event_Id = evtHeader.GetId();
      Run_Id = evtHeader.GetRunNumber();
      const det::Detector& detector = det::Detector::GetInstance();
      //detector.Update(Event_Timestamp,Run_Id);

      const det::MagneticFieldTracker& tracker = detector.GetMagneticFieldTracker();

      //SimEvent
        using namespace evt::sim;
        const SimEvent& simEvent = event.GetSimEvent();
        if (event.IsSimulation())
        {
        cout << "FlowQA::Calculate_QA //SimEvent" << endl;
	//Event
	    Sim_ImpactParameter = simEvent.GetPrimaryInteraction().GetImpactParameter();
	    Sim_ReactionPlaneAngle = simEvent.GetPrimaryInteraction().GetEventPlaneAngle();
        //Tracks
            Sim_nTracks = 0;
            int counter = 0;
            if (simEvent.HasMainVertex())
            {
                for (std::list<evt::sim::VertexTrack>::const_iterator simTrackIter = simEvent.Begin<evt::sim::VertexTrack>(), simTrackEnd = simEvent.End<evt::sim::VertexTrack>(); simTrackIter != simTrackEnd; ++simTrackIter)
                {

                    cout << "FlowQA::Calculate_QA //SimEvent counter = " << counter << endl;
                    counter++;
//                     const evt::sim::VertexTrack& simTrack = simEvent.Get(*vtxTrackIter);
                    const evt::sim::VertexTrack& simTrack = *simTrackIter;
                    const Vector& momentum = simTrack.GetMomentum();
                    Sim_track_pT[counter] = momentum.GetRho();
                    Sim_track_phi[counter] = momentum.GetPhi();
                    Sim_track_eta[counter] = 0.5*TMath::Log((momentum.GetR()+momentum.GetZ())/(momentum.GetR()-momentum.GetZ()));
//                     Sim_track_mass[counter] = ???;
                    Sim_track_charge[counter] = simTrack.GetCharge();
                    Sim_track_pdg_id[counter] = simTrack.GetParticleId();
                    Sim_track_IsInTPC_Total[counter] = (simTrack.IsInTPC(det::TPCConst::eVTPC1)+simTrack.IsInTPC(det::TPCConst::eVTPC2)+simTrack.IsInTPC(det::TPCConst::eMTPCL)+simTrack.IsInTPC(det::TPCConst::eMTPCR) > 0);
                    Sim_track_IsInTPC_TPCV1[counter] = simTrack.IsInTPC(det::TPCConst::eVTPC1);
                    Sim_track_IsInTPC_TPCV2[counter] = simTrack.IsInTPC(det::TPCConst::eVTPC2);
                    Sim_track_IsInTPC_TPCVmain[counter] = (simTrack.IsInTPC(det::TPCConst::eMTPCL)+simTrack.IsInTPC(det::TPCConst::eMTPCR) > 0);
                    
                    Sim_track_Number_of_Hits[counter] = simTrack.GetNumberOfHits();
                    Sim_track_IsSpectator[counter] = simTrack.IsSpectator();
                    Sim_nTracks++;
//         bool Sim_track_IsInTPC[nMaxTracks];
                }
            }

        //PSD
            evt::sim::PSD& simPSD = event.GetSimEvent().GetPSD();
            for (int i=0;i<simPSD.GetNModules();i++)
            {
                evt::sim::PSDModule simPSDModule = simPSD.GetModule(i+1);
                Sim_PSD_module_Temperature[i] = simPSDModule.GetTemperature();
                for (int j=0;j<simPSDModule.GetNSections();j++)
                {
//                     cout << "mod " << i << " sec " << j << endl;
                    evt::sim::PSDSection simPSDSecection = simPSDModule.GetSection(j+1);
                    Sim_PSD_Energy += simPSDSecection.GetEnergy();
                    Sim_PSD_module_Energy[i] += simPSDSecection.GetEnergy();
                    Sim_PSD_section_Energy[i][j] = simPSDSecection.GetEnergy();
                }
            }
        }
      //RecEvent
        const evt::RecEvent& recEvent = event.GetRecEvent();
        if (!event.IsSimulation())
        {
            cout << "FlowQA::Calculate_QA //RecEvent" << endl;
            //Main vertex
            const evt::rec::Vertex& mainVertex = recEvent.GetMainVertex();
            const Point& pm = mainVertex.GetPosition();
            Main_Vertex_X = pm.GetX();
            Main_Vertex_Y = pm.GetY();
            Main_Vertex_Z = pm.GetZ();
            if (event.GetRecEvent().HasPrimaryVertex(rec::VertexConst::ePrimaryFitZ)){
                const evt::rec::Vertex& pVtx = recEvent.GetPrimaryVertex(rec::VertexConst::ePrimaryFitZ);
                const Point& pntVtx = pVtx.GetPosition();
                pVtx_Z = pntVtx.GetZ();
            } else{
                pVtx_Z = -999;
            } 
           //TPC
                float px;
                float py;
                float pz;
                float p;
                float pT;
                float eta;
                float phi;
                int charge;
                Point endPos, endPosFirst, endPosLast, endTestPos, endTestPosFirst, endTestPosLast, endPSDPos;
                Vector endP, endPFirst, endPLast, endTestP, endTestPFirst, endTestPLast, endPSDP;
				const float PSD_Z = 17000.; 
		
	   //Main vertex tracks
		int track_id = 0;
                TPC_main_vtx_nTracks = 0;
                for (evt::rec::VertexTrackIndexIterator vtxTrackIter = mainVertex.DaughterTracksBegin(); vtxTrackIter != mainVertex.DaughterTracksEnd(); ++vtxTrackIter)
                {
                    TPC_Multiplicity_all++;
                    const evt::rec::VertexTrack& vtxTrack = recEvent.Get(*vtxTrackIter);
		    if (! (vtxTrack.GetStatus() == 0 && vtxTrack.HasTrack())) continue;
                    TPC_Multiplicity++;
		    //main vtx track (momentum of track from main vertex (for pA reconstruciton this is BPD vertex, nominal) at the vertex point)
		    px=vtxTrack.GetMomentum().GetX();
                    py=vtxTrack.GetMomentum().GetY();
                    pz=vtxTrack.GetMomentum().GetZ();
                    p=TMath::Sqrt(px*px+py*py+pz*pz);
                    pT=TMath::Sqrt(px*px+py*py);
                    eta = 0.5*TMath::Log((p+pz)/(p-pz));
                    phi = TMath::ATan2(py,px);
                    if (TMath::IsNaN(eta)) { cout << "ERROR fTreeQA: eta is NaN!" << endl; continue; }
                    if (TMath::IsNaN(phi)) { cout << "ERROR fTreeQA: phi is NaN!" << endl; continue; }
                    TPC_main_vtx_pT[track_id] = pT;
                    TPC_main_vtx_px[track_id] = px;
                    TPC_main_vtx_py[track_id] = py;
                    TPC_main_vtx_pz[track_id] = pz;
                    TPC_main_vtx_eta[track_id] = eta;
                    TPC_main_vtx_phi[track_id] = phi;
                    //get tracks in TPC associated with the main vertex
		    const evt::rec::Track& track = recEvent.Get(vtxTrack.GetTrackIndex());
                    using namespace TrackConst;
		    //track momentum at first point in TPC (track from main vertex)
                    px = track.GetMomentum().GetX();
                    py = track.GetMomentum().GetY();
                    pz = track.GetMomentum().GetZ();
                    p=TMath::Sqrt(px*px+py*py+pz*pz);
                    pT=TMath::Sqrt(px*px+py*py);
                    eta = 0.5*TMath::Log((p+pz)/(p-pz));
                    phi = TMath::ATan2(py,px);
                    if (TMath::IsNaN(eta)) { cout << "ERROR fTreeQA: eta is NaN!" << endl; continue; }
                    if (TMath::IsNaN(phi)) { cout << "ERROR fTreeQA: phi is NaN!" << endl; continue; }
                    TPC_main_vtx_track_pT[track_id] = pT;
                    TPC_main_vtx_track_px[track_id] = px;
                    TPC_main_vtx_track_py[track_id] = py;
                    TPC_main_vtx_track_pz[track_id] = pz;
                    TPC_main_vtx_track_eta[track_id] = eta;
                    TPC_main_vtx_track_phi[track_id] = phi;
                    TPC_main_vtx_track_point_x[track_id] = track.GetMomentumPoint().GetX();
                    TPC_main_vtx_track_point_y[track_id] = track.GetMomentumPoint().GetY();
                    TPC_main_vtx_track_point_z[track_id] = track.GetMomentumPoint().GetZ();
                    TPC_main_vtx_track_first_point_x[track_id] = track.GetFirstPointOnTrack().GetX();
                    TPC_main_vtx_track_first_point_y[track_id] = track.GetFirstPointOnTrack().GetY();
                    TPC_main_vtx_track_first_point_z[track_id] = track.GetFirstPointOnTrack().GetZ();
                    TPC_main_vtx_track_last_point_x[track_id] = track.GetLastPointOnTrack().GetX();
                    TPC_main_vtx_track_last_point_y[track_id] = track.GetLastPointOnTrack().GetY();
                    TPC_main_vtx_track_last_point_z[track_id] = track.GetLastPointOnTrack().GetZ();
                    
					if (track.GetMomentumPoint().GetZ()) tracker.TrackToZ(track.GetMomentumPoint().GetZ(), track.GetCharge(), vtxTrack.GetImpactPoint(), vtxTrack.GetMomentum(), endTestPos, endTestP);
                    if (track.GetFirstPointOnTrack().GetZ()) tracker.TrackToZ(track.GetFirstPointOnTrack().GetZ(), track.GetCharge(), vtxTrack.GetImpactPoint(), vtxTrack.GetMomentum(), endTestPosFirst, endTestPFirst);
                    if (track.GetLastPointOnTrack().GetZ()) tracker.TrackToZ(track.GetLastPointOnTrack().GetZ(), track.GetCharge(), vtxTrack.GetImpactPoint(), vtxTrack.GetMomentum(), endTestPosLast, endTestPLast);
                    //if (pVtx_Z == -999) pVtx_Z = vtxTrack.GetImpactPoint().GetZ();
					if (pVtx_Z != -999){
                        if (track.GetMomentumPoint().GetZ()) tracker.TrackToZ(pVtx_Z, vtxTrack.GetCharge(), track.GetMomentumPoint(), track.GetMomentum(), endPos, endP);
                        if (track.GetFirstPointOnTrack().GetZ()) tracker.TrackToZ(pVtx_Z, vtxTrack.GetCharge(), track.GetFirstPointOnTrack(), track.GetMomentum(), endPosFirst, endPFirst);
                        if (track.GetLastPointOnTrack().GetZ()) tracker.TrackToZ(pVtx_Z, vtxTrack.GetCharge(), track.GetLastPointOnTrack(), track.GetMomentum(), endPosLast, endPLast);
                    }
					
					tracker.TrackToZ(PSD_Z, track.GetCharge(), vtxTrack.GetImpactPoint(), vtxTrack.GetMomentum(), endPSDPos, endPSDP);

                    TPC_main_vtx_impact_point_X[track_id] = vtxTrack.GetImpactPoint().GetX();
                    TPC_main_vtx_impact_point_Y[track_id] = vtxTrack.GetImpactPoint().GetY();
                    TPC_main_vtx_impact_point_Z[track_id] = vtxTrack.GetImpactPoint().GetZ();
                    TPC_main_vtx_ext_track_px[track_id] = endP.GetX();
                    TPC_main_vtx_ext_track_py[track_id] = endP.GetY();
                    TPC_main_vtx_ext_track_pz[track_id] = endP.GetZ();
                    TPC_main_vtx_ext_first_track_px[track_id] = endPFirst.GetX();
                    TPC_main_vtx_ext_first_track_py[track_id] = endPFirst.GetY();
                    TPC_main_vtx_ext_first_track_pz[track_id] = endPFirst.GetZ();
                    TPC_main_vtx_ext_last_track_px[track_id] = endPLast.GetX();
                    TPC_main_vtx_ext_last_track_py[track_id] = endPLast.GetY();
                    TPC_main_vtx_ext_last_track_pz[track_id] = endPLast.GetZ();
                    TPC_main_vtx_ext_point_x[track_id] = endPos.GetX();
                    TPC_main_vtx_ext_point_y[track_id] = endPos.GetY();
                    TPC_main_vtx_ext_point_z[track_id] = endPos.GetZ();
                    TPC_main_vtx_ext_first_point_x[track_id] = endPosFirst.GetX();
                    TPC_main_vtx_ext_first_point_y[track_id] = endPosFirst.GetY();
                    TPC_main_vtx_ext_first_point_z[track_id] = endPosFirst.GetZ();
                    TPC_main_vtx_ext_last_point_x[track_id] = endPosLast.GetX();
                    TPC_main_vtx_ext_last_point_y[track_id] = endPosLast.GetY();
                    TPC_main_vtx_ext_last_point_z[track_id] = endPosLast.GetZ();

					TPC_main_test_ext_track_px[track_id] = endTestP.GetX();
                    TPC_main_test_ext_track_py[track_id] = endTestP.GetY();
                    TPC_main_test_ext_track_pz[track_id] = endTestP.GetZ();
                    TPC_main_test_ext_first_track_px[track_id] = endTestPFirst.GetX();
                    TPC_main_test_ext_first_track_py[track_id] = endTestPFirst.GetY();
                    TPC_main_test_ext_first_track_pz[track_id] = endTestPFirst.GetZ();
                    TPC_main_test_ext_last_track_px[track_id] = endTestPLast.GetX();
                    TPC_main_test_ext_last_track_py[track_id] = endTestPLast.GetY();
                    TPC_main_test_ext_last_track_pz[track_id] = endTestPLast.GetZ();
                    TPC_main_test_ext_point_x[track_id] = endTestPos.GetX();
                    TPC_main_test_ext_point_y[track_id] = endTestPos.GetY();
                    TPC_main_test_ext_point_z[track_id] = endTestPos.GetZ();
                    TPC_main_test_ext_first_point_x[track_id] = endTestPosFirst.GetX();
                    TPC_main_test_ext_first_point_y[track_id] = endTestPosFirst.GetY();
                    TPC_main_test_ext_first_point_z[track_id] = endTestPosFirst.GetZ();
                    TPC_main_test_ext_last_point_x[track_id] = endTestPosLast.GetX();
                    TPC_main_test_ext_last_point_y[track_id] = endTestPosLast.GetY();
                    TPC_main_test_ext_last_point_z[track_id] = endTestPosLast.GetZ();

					TPC_main_PSD_ext_track_px[track_id] = endPSDP.GetX();
                    TPC_main_PSD_ext_track_py[track_id] = endPSDP.GetY();
                    TPC_main_PSD_ext_track_pz[track_id] = endPSDP.GetZ();
					TPC_main_PSD_ext_point_x[track_id] = endPSDPos.GetX();
                    TPC_main_PSD_ext_point_y[track_id] = endPSDPos.GetY();
                    TPC_main_PSD_ext_point_z[track_id] = endPSDPos.GetZ();


                    //track TPC parameters
		            TPC_main_vtx_track_charge[track_id] = track.GetCharge();
                    TPC_main_vtx_track_nClusters_Total[track_id] = track.GetNumberOfClusters(eAll);
                    TPC_main_vtx_track_nClusters_TPCV1[track_id] = track.GetNumberOfClusters(eVTPC1);
                    TPC_main_vtx_track_nClusters_TPCV2[track_id] = track.GetNumberOfClusters(eVTPC2);
                    TPC_main_vtx_track_nClusters_TPCVmain[track_id] = track.GetNumberOfClusters(eMTPC);
                    TPC_main_vtx_track_nClusters_TPCVgap[track_id] = track.GetNumberOfClusters(eGTPC);
                    TPC_main_vtx_track_nClustersPotential_Total[track_id] = track.GetPotentialNumberOfClusters(eAll);
                    TPC_main_vtx_track_nClustersPotential_TPCV1[track_id] = track.GetPotentialNumberOfClusters(eVTPC1);
                    TPC_main_vtx_track_nClustersPotential_TPCV2[track_id] = track.GetPotentialNumberOfClusters(eVTPC2);
                    TPC_main_vtx_track_nClustersPotential_TPCVmain[track_id] = track.GetPotentialNumberOfClusters(eMTPC);
                    TPC_main_vtx_track_nClustersPotential_TPCVgap[track_id] = track.GetPotentialNumberOfClusters(eGTPC);
                    TPC_main_vtx_track_nClustersFit_Total[track_id] = track.GetNumberOfFitClusters(eAll);
                    TPC_main_vtx_track_nClustersFit_TPCV1[track_id] = track.GetNumberOfFitClusters(eVTPC1);
                    TPC_main_vtx_track_nClustersFit_TPCV2[track_id] = track.GetNumberOfFitClusters(eVTPC2);
                    TPC_main_vtx_track_nClustersFit_TPCVmain[track_id] = track.GetNumberOfFitClusters(eMTPC);
                    TPC_main_vtx_track_nClustersFit_TPCVgap[track_id] = track.GetNumberOfFitClusters(eGTPC);
                    TPC_main_vtx_track_nClustersdEdX_Total[track_id] = track.GetNumberOfdEdXClusters(eAll);
                    TPC_main_vtx_track_nClustersdEdX_TPCV1[track_id] = track.GetNumberOfdEdXClusters(eVTPC1);
                    TPC_main_vtx_track_nClustersdEdX_TPCV2[track_id] = track.GetNumberOfdEdXClusters(eVTPC2);
                    TPC_main_vtx_track_nClustersdEdX_TPCVmain[track_id] = track.GetNumberOfdEdXClusters(eMTPC);
                    TPC_main_vtx_track_nClustersdEdX_TPCVgap[track_id] = track.GetNumberOfdEdXClusters(eGTPC);
                    TPC_main_vtx_track_EnergyClusters_Total[track_id] = track.GetEnergyDeposit(eAll);
                    TPC_main_vtx_track_EnergyClusters_TPCV1[track_id] = track.GetEnergyDeposit(eVTPC1);
                    TPC_main_vtx_track_EnergyClusters_TPCV2[track_id] = track.GetEnergyDeposit(eVTPC2);
                    TPC_main_vtx_track_EnergyClusters_TPCVmain[track_id] = track.GetEnergyDeposit(eMTPC);
                    TPC_main_vtx_track_EnergyClusters_TPCVgap[track_id] = track.GetEnergyDeposit(eGTPC);
                    TPC_main_vtx_track_chi2[track_id] = vtxTrack.GetChi2();
                    TPC_main_vtx_track_ndf[track_id] = vtxTrack.GetNdf();
                    track_id++;
                    TPC_main_vtx_nTracks++;
                    if (track_id > nMaxTracks)
                    {
                        for (int k=0;k<nMaxTracks;k++)
                        {
                            TPC_main_vtx_track_pT[k] = undefinedValue;
                        }
                    }
                }
            
            Has_Primary_Vertex = false;    
	    if (event.GetRecEvent().HasPrimaryVertex(rec::VertexConst::ePrimaryFitZ))
	    {
		//Primary vertex
		Has_Primary_Vertex = true;
		const evt::rec::Vertex& primaryVertex = recEvent.GetPrimaryVertex(rec::VertexConst::ePrimaryFitZ);
		const Point& pm = primaryVertex.GetPosition();
		Primary_Vertex_X = pm.GetX();
		Primary_Vertex_Y = pm.GetY();
		Primary_Vertex_Z = pm.GetZ();
		Primary_Vertex_Q = (int) primaryVertex.GetFitQuality();
		//Primary vertex tracks
		int track_id = 0;
                TPC_primary_vtx_nTracks = 0;
                for (evt::rec::VertexTrackIndexIterator vtxTrackIter = primaryVertex.DaughterTracksBegin(); vtxTrackIter != primaryVertex.DaughterTracksEnd(); ++vtxTrackIter)
                {
                    const evt::rec::VertexTrack& vtxTrack = recEvent.Get(*vtxTrackIter);
		    if (! (vtxTrack.GetStatus() == 0 && vtxTrack.HasTrack())) continue;
		    //priamry vtx track (momentum of track from primary vertex (fitted vertex position) at the vertex point)
		    px=vtxTrack.GetMomentum().GetX();
                    py=vtxTrack.GetMomentum().GetY();
                    pz=vtxTrack.GetMomentum().GetZ();
                    p=TMath::Sqrt(px*px+py*py+pz*pz);
                    pT=TMath::Sqrt(px*px+py*py);
                    eta = 0.5*TMath::Log((p+pz)/(p-pz));
                    phi = TMath::ATan2(py,px);
                    if (TMath::IsNaN(eta)) { cout << "ERROR fTreeQA: eta is NaN!" << endl; continue; }
                    if (TMath::IsNaN(phi)) { cout << "ERROR fTreeQA: phi is NaN!" << endl; continue; }
                    TPC_primary_vtx_pT[track_id] = pT;
                    TPC_primary_vtx_eta[track_id] = eta;
                    TPC_primary_vtx_phi[track_id] = phi;
                    //get tracks in TPC associated with the primary vertex
		    const evt::rec::Track& track = recEvent.Get(vtxTrack.GetTrackIndex());
                    using namespace TrackConst;
		    //track momentum at first point in TPC (track from primary vertex)
                    px = track.GetMomentum().GetX();
                    py = track.GetMomentum().GetY();
                    pz = track.GetMomentum().GetZ();
                    p=TMath::Sqrt(px*px+py*py+pz*pz);
                    pT=TMath::Sqrt(px*px+py*py);
                    eta = 0.5*TMath::Log((p+pz)/(p-pz));
                    phi = TMath::ATan2(py,px);
                    if (TMath::IsNaN(eta)) { cout << "ERROR fTreeQA: eta is NaN!" << endl; continue; }
                    if (TMath::IsNaN(phi)) { cout << "ERROR fTreeQA: phi is NaN!" << endl; continue; }
                    TPC_primary_vtx_track_pT[track_id] = pT;
                    TPC_primary_vtx_track_eta[track_id] = eta;
                    TPC_primary_vtx_track_phi[track_id] = phi;
		    TPC_primary_vtx_impact_point_X[track_id] = vtxTrack.GetImpactPoint().GetX();
		    TPC_primary_vtx_impact_point_Y[track_id] = vtxTrack.GetImpactPoint().GetY();
		    TPC_primary_vtx_impact_point_Z[track_id] = vtxTrack.GetImpactPoint().GetZ();
                    //track TPC parameters
		    TPC_primary_vtx_track_charge[track_id] = track.GetCharge();
                    TPC_primary_vtx_track_nClusters_Total[track_id] = track.GetNumberOfClusters(eAll);
                    TPC_primary_vtx_track_nClusters_TPCV1[track_id] = track.GetNumberOfClusters(eVTPC1);
                    TPC_primary_vtx_track_nClusters_TPCV2[track_id] = track.GetNumberOfClusters(eVTPC2);
                    TPC_primary_vtx_track_nClusters_TPCVmain[track_id] = track.GetNumberOfClusters(eMTPC);
                    TPC_primary_vtx_track_nClusters_TPCVgap[track_id] = track.GetNumberOfClusters(eGTPC);
                    TPC_primary_vtx_track_nClustersPotential_Total[track_id] = track.GetPotentialNumberOfClusters(eAll);
                    TPC_primary_vtx_track_nClustersPotential_TPCV1[track_id] = track.GetPotentialNumberOfClusters(eVTPC1);
                    TPC_primary_vtx_track_nClustersPotential_TPCV2[track_id] = track.GetPotentialNumberOfClusters(eVTPC2);
                    TPC_primary_vtx_track_nClustersPotential_TPCVmain[track_id] = track.GetPotentialNumberOfClusters(eMTPC);
                    TPC_primary_vtx_track_nClustersPotential_TPCVgap[track_id] = track.GetPotentialNumberOfClusters(eGTPC);
                    TPC_primary_vtx_track_nClustersFit_Total[track_id] = track.GetNumberOfFitClusters(eAll);
                    TPC_primary_vtx_track_nClustersFit_TPCV1[track_id] = track.GetNumberOfFitClusters(eVTPC1);
                    TPC_primary_vtx_track_nClustersFit_TPCV2[track_id] = track.GetNumberOfFitClusters(eVTPC2);
                    TPC_primary_vtx_track_nClustersFit_TPCVmain[track_id] = track.GetNumberOfFitClusters(eMTPC);
                    TPC_primary_vtx_track_nClustersFit_TPCVgap[track_id] = track.GetNumberOfFitClusters(eGTPC);
                    TPC_primary_vtx_track_nClustersdEdX_Total[track_id] = track.GetNumberOfdEdXClusters(eAll);
                    TPC_primary_vtx_track_nClustersdEdX_TPCV1[track_id] = track.GetNumberOfdEdXClusters(eVTPC1);
                    TPC_primary_vtx_track_nClustersdEdX_TPCV2[track_id] = track.GetNumberOfdEdXClusters(eVTPC2);
                    TPC_primary_vtx_track_nClustersdEdX_TPCVmain[track_id] = track.GetNumberOfdEdXClusters(eMTPC);
                    TPC_primary_vtx_track_nClustersdEdX_TPCVgap[track_id] = track.GetNumberOfdEdXClusters(eGTPC);
                    TPC_primary_vtx_track_EnergyClusters_Total[track_id] = track.GetEnergyDeposit(eAll);
                    TPC_primary_vtx_track_EnergyClusters_TPCV1[track_id] = track.GetEnergyDeposit(eVTPC1);
                    TPC_primary_vtx_track_EnergyClusters_TPCV2[track_id] = track.GetEnergyDeposit(eVTPC2);
                    TPC_primary_vtx_track_EnergyClusters_TPCVmain[track_id] = track.GetEnergyDeposit(eMTPC);
                    TPC_primary_vtx_track_EnergyClusters_TPCVgap[track_id] = track.GetEnergyDeposit(eGTPC);
                    TPC_primary_vtx_track_chi2[track_id] = vtxTrack.GetChi2();
                    TPC_primary_vtx_track_ndf[track_id] = vtxTrack.GetNdf();
                    track_id++;
                    TPC_primary_vtx_nTracks++;
                    if (track_id > nMaxTracks)
                    {
                        std::cout << "Too many tracks!" << std::endl;
                        for (int k=0;k<nMaxTracks;k++)
                        {
                            TPC_primary_vtx_track_pT[k] = undefinedValue;
                        }
                    }
                }
	    }
            //PSD
            const evt::rec::PSD& recPSDEvent = recEvent.GetPSD();
            for (int i=0; i<(int)recPSDEvent.GetNModules(); i++)
            {
                evt::rec::PSDModule fPSDmodule = recPSDEvent.GetModule(i+1);
                PSD_module_number_of_sections[i] = fPSDmodule.GetNSections();
                PSD_module_Energy_default[i] = fPSDmodule.GetEnergy();
                PSD_module_X[i] = detector.GetPSD().GetPSDModule(det::PSDConst::GetType(det::PSDConst::GetIdByNumber(i+1)))->GetModulePositionX(det::PSDConst::GetIdByNumber(i+1));
                PSD_module_Y[i] = detector.GetPSD().GetPSDModule(det::PSDConst::GetType(det::PSDConst::GetIdByNumber(i+1)))->GetModulePositionY(det::PSDConst::GetIdByNumber(i+1));
                for (int j=0; j<(int)fPSDmodule.GetNSections(); j++)
                {
                    evt::rec::PSDSection fPSDsection = fPSDmodule.GetSection(j+1);
                    PSD_section_Energy[i][j] = fPSDsection.GetEnergy();
                    if (i!=44) PSD_section_slice_Energy[j] += PSD_section_Energy[i][j];
                    PSD_module_Energy[i] += fPSDsection.GetEnergy();
                    if (PSD_section_Energy[i][j] > 0) {PSD_section_Number++;PSD_section_Number_array[i][j]++;}
                }
                if (PSD_module_Energy[i] > 0) {PSD_module_Number++;PSD_module_Number_array[i]++;}
                if (i<16 || i==44) PSD_1_Energy += PSD_module_Energy[i];
                if (i>=16 && i<28) PSD_2_Energy += PSD_module_Energy[i];
                if (i>=28 && i<44) PSD_3_Energy += PSD_module_Energy[i];
                PSD_Energy += PSD_module_Energy[i];
            }
        }
      //Triggers (old style variables)
      T1 = trigBit[eTrig1];
      T2 = trigBit[eTrig2];
      T3 = trigBit[eTrig3];
      T4 = trigBit[eTrig4];
      
      
      //BEGIN TEST
      //Triggers
      
      if (!event.IsSimulation())
      {
            evt::rec::BPDPlane fBPDPlaneReco[nBPD][nBPDcomponents];
            fBPDPlaneReco[0][0] = event.GetRecEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD1, det::BPDConst::eX);
            fBPDPlaneReco[0][1] = event.GetRecEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD1, det::BPDConst::eY);
            fBPDPlaneReco[1][0] = event.GetRecEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD2, det::BPDConst::eX);
            fBPDPlaneReco[1][1] = event.GetRecEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD2, det::BPDConst::eY);
            fBPDPlaneReco[2][0] = event.GetRecEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD3, det::BPDConst::eX);
            fBPDPlaneReco[2][1] = event.GetRecEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD3, det::BPDConst::eY);
      
        //       evt::raw::BPDPlane fBPDPlaneRaw[nBPD][nBPDcomponents];
        //       fBPDPlaneRaw[0][0] = event.GetRawEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD1, det::BPDConst::eX);
        //       fBPDPlaneRaw[0][1] = event.GetRawEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD1, det::BPDConst::eY);
        //       fBPDPlaneRaw[1][0] = event.GetRawEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD2, det::BPDConst::eX);
        //       fBPDPlaneRaw[1][1] = event.GetRawEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD2, det::BPDConst::eY);
        //       fBPDPlaneRaw[2][0] = event.GetRawEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD3, det::BPDConst::eX);
        //       fBPDPlaneRaw[2][1] = event.GetRawEvent().GetBeam().GetBPDPlane(det::BPDConst::eBPD3, det::BPDConst::eY);
            
            for (int i=0;i<nBPD;i++)
            {
                for (int j=0;j<nBPDcomponents;j++)
                {
        // 	      BPD_Status[0][i][j] = fBPDPlaneRaw[i][j].GetStatus();
        // 	      BPD_Position[0][i][j] = fBPDPlaneRaw[i][j].GetPosition();
        // 	      BPD_PositionError[0][i][j] = fBPDPlaneRaw[i][j].GetPositionError();
        // 	      BPD_Z[0][i][j] = fBPDPlaneRaw[i][j].GetZ();
        // 	      BPD_RMS[0][i][j] = fBPDPlaneRaw[i][j].GetRMS();
        // 	      BPD_Maximum[0][i][j] = fBPDPlaneRaw[i][j].GetMaximum();
        // 	      BPD_Charge[0][i][j] = fBPDPlaneRaw[i][j].GetCharge();
        // 	      BPD_SumOfAll[0][i][j] = fBPDPlaneRaw[i][j].GetSumOfAll();
                    
                    BPD_Status[1][i][j] = fBPDPlaneReco[i][j].GetStatus();
                    BPD_Position[1][i][j] = fBPDPlaneReco[i][j].GetPosition();
                    BPD_PositionError[1][i][j] = fBPDPlaneReco[i][j].GetPositionError();
                    BPD_Z[1][i][j] = fBPDPlaneReco[i][j].GetZ();
                    BPD_RMS[1][i][j] = fBPDPlaneReco[i][j].GetRMS();
                    BPD_Maximum[1][i][j] = fBPDPlaneReco[i][j].GetMaximum();
                    BPD_Charge[1][i][j] = fBPDPlaneReco[i][j].GetCharge();
                    BPD_SumOfAll[1][i][j] = fBPDPlaneReco[i][j].GetSumOfAll();
                }
            }
//             evt::rec::Trigger fTriggerRec = event.GetRecEvent().GetBeam().GetTrigger();
      }
      evt::rec::Trigger fTriggerRec = event.GetRecEvent().GetBeam().GetTrigger();
      evt::raw::Trigger fTriggerRaw = event.GetRawEvent().GetBeam().GetTrigger();
      
//       const evt::raw::Trigger& trigger = event.GetRawEvent().GetBeam().GetTrigger();
//       const unsigned int V1Padc = trigger.GetADC( TriggerConst::eV1p );
//       fQAHistograms1D[eV1PADC].Fill( double(V1Padc) );
      
      triggersADC[0][0] = fTriggerRaw.GetADC(det::TriggerConst::eS1);
      triggersADC[0][1] = fTriggerRaw.GetADC(det::TriggerConst::eS2);
      triggersADC[0][2] = fTriggerRaw.GetADC(det::TriggerConst::eS3);
      triggersADC[0][3] = fTriggerRaw.GetADC(det::TriggerConst::eV1);
      triggersADC[0][4] = fTriggerRaw.GetADC(det::TriggerConst::eV1p);
      triggersADC[0][5] = fTriggerRaw.GetADC(det::TriggerConst::ePSD);

//       const evt::raw::Trigger& fTriggerRaw1 = event.GetRawEvent().GetBeam().GetTrigger();
//       cout << "DEBUG: without &: " << fTriggerRaw1.GetADC(det::TriggerConst::eV1p) << "; with &: " << fTriggerRaw.GetADC(det::TriggerConst::eV1p) 
//       << " triggersADC[0][4] " << triggersADC[0][4] << endl;

//       triggersADC[0][0] = float(fTriggerRaw.GetADC(det::TriggerConst::eS1));
//       triggersADC[0][1] = float(fTriggerRaw.GetADC(det::TriggerConst::eS2));
//       triggersADC[0][2] = float(fTriggerRaw.GetADC(det::TriggerConst::eS3));
//       triggersADC[0][3] = float(fTriggerRaw.GetADC(det::TriggerConst::eV1));
//       triggersADC[0][4] = float(fTriggerRaw.GetADC(det::TriggerConst::eV1p));
//       triggersADC[0][5] = float(fTriggerRaw.GetADC(det::TriggerConst::ePSD));

//       triggersADC[1][0] = fTriggerRec.GetCharge(det::TriggerConst::eS1);
//       triggersADC[1][1] = fTriggerRec.GetCharge(det::TriggerConst::eS2);
//       triggersADC[1][2] = fTriggerRec.GetCharge(det::TriggerConst::eS3);
//       triggersADC[1][3] = fTriggerRec.GetCharge(det::TriggerConst::eV1);
//       triggersADC[1][4] = fTriggerRec.GetCharge(det::TriggerConst::eV1p);
//       triggersADC[1][5] = fTriggerRec.GetCharge(det::TriggerConst::ePSD);
      
      isTriggers_Simple[0][0] = fTriggerRaw.IsTrigger(det::TriggerConst::eS1,det::TriggerConst::ePrescaled);
      isTriggers_Simple[0][1] = fTriggerRaw.IsTrigger(det::TriggerConst::eS2,det::TriggerConst::ePrescaled);
      isTriggers_Simple[0][2] = fTriggerRaw.IsTrigger(det::TriggerConst::eS3,det::TriggerConst::ePrescaled);
      isTriggers_Simple[0][3] = fTriggerRaw.IsTrigger(det::TriggerConst::eV1,det::TriggerConst::ePrescaled);
      isTriggers_Simple[0][4] = fTriggerRaw.IsTrigger(det::TriggerConst::eV1p,det::TriggerConst::ePrescaled);
      isTriggers_Simple[0][5] = fTriggerRaw.IsTrigger(det::TriggerConst::ePSD,det::TriggerConst::ePrescaled);
      
      isTriggers_Combined[0][0] = fTriggerRaw.IsTrigger(det::TriggerConst::eT1,det::TriggerConst::ePrescaled);
      isTriggers_Combined[0][1] = fTriggerRaw.IsTrigger(det::TriggerConst::eT2,det::TriggerConst::ePrescaled);
      isTriggers_Combined[0][2] = fTriggerRaw.IsTrigger(det::TriggerConst::eT3,det::TriggerConst::ePrescaled);
      isTriggers_Combined[0][3] = fTriggerRaw.IsTrigger(det::TriggerConst::eT4,det::TriggerConst::ePrescaled);

      if (!event.IsSimulation())
      {
	  Beam_Momentum[1][0] = event.GetRecEvent().GetBeam().GetMomentum().GetX();
	  Beam_Momentum[1][1] = event.GetRecEvent().GetBeam().GetMomentum().GetY();
	  Beam_Momentum[1][2] = event.GetRecEvent().GetBeam().GetMomentum().GetZ();
	  Beam_Status[1] = event.GetRecEvent().GetBeam().GetStatus();
      }
      //END TEST
      
      if (!event.IsSimulation())
      {
          if (fTriggerRaw.HasTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1_1))
          {
                WFA_NumberOfSignalHits[0] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::eS1_1);
                for (int j=0;j<nMaxWFAsignals;j++)
                {
                    WFA_TimeStructure[0][j] = undefinedValue;
                }
                for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1_1).size();j++)
                {
                    WFA_TimeStructure[0][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1_1).at(j);
                }
	  }
          if (fTriggerRaw.HasTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eT4))
          {
                WFA_NumberOfSignalHits[1] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::eT4);
                for (int j=0;j<nMaxWFAsignals;j++)
                {
                    WFA_TimeStructure[1][j] = undefinedValue;
                }
                for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eT4).size();j++)
                {
                    WFA_TimeStructure[1][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eT4).at(j);
                }
	  }
      }
      //timestructure
//             det::TriggerConst::EId TriggerConsts[6];
//             TriggerConsts[0] = det::TriggerConst::eS1;
//             TriggerConsts[1] = det::TriggerConst::eS2;
//             TriggerConsts[2] = det::TriggerConst::eS3;
//             TriggerConsts[3] = det::TriggerConst::eV1;
//             TriggerConsts[4] = det::TriggerConst::eV1p;
//             TriggerConsts[5] = det::TriggerConst::ePSD;
            
//             for (int i=0;i<6;i++)
//             {
//                 cout << "TriggerRaw[" << i << "] has time structure = " << fTriggerRaw.HasTimeStructure(det::TimeStructureConst::eWFA,TriggerConsts[i]) << endl;
//                 if (fTriggerRaw.HasTimeStructure(det::TimeStructureConst::eWFA,TriggerConsts[i]))
//                 {
//                     WFA_NumberOfSignalHits[i] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,TriggerConsts[i]);
//                     cout << fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,TriggerConsts[i]) << endl;
//                 }
//                 cout << "TriggerRec[" << i << "] has time structure = " << fTriggerRec.HasTimeStructure(det::TimeStructureConst::eWFA,TriggerConsts[i]) << endl;
//                 if (fTriggerRec.HasTimeStructure(det::TimeStructureConst::eWFA,TriggerConsts[i]))
//                 {
//                     WFA_NumberOfSignalHits[i] = fTriggerRec.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,TriggerConsts[i]);
//                     cout << fTriggerRec.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,TriggerConsts[i]) << endl;
//                 }
//             }
// 	}
//       WFA_NumberOfSignalHits[0] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::eS1);
//       WFA_NumberOfSignalHits[1] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::eS2);
//       WFA_NumberOfSignalHits[2] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::eS3);
//       WFA_NumberOfSignalHits[3] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::eV1);
//       WFA_NumberOfSignalHits[4] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::eV1p);
//       WFA_NumberOfSignalHits[5] = fTriggerRaw.GetNumberOfSignalHits(det::TimeStructureConst::eWFA,det::TriggerConst::ePSD);
      
//       //S1
//       cout << "test 0" << endl;
//       cout << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(0) << endl;
//       cout << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(1) << endl;
//       cout << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(2) << endl;
//       cout << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(3) << endl;
//       cout << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(4) << endl;
//       cout << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(5) << endl;
//       cout << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).size() << endl;
//       cout << "test 1" << endl;
//       if (fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).size() > nMaxWFAsignals)
//       {
//           cout << "test 1 1" << endl;
// 	  for (int j=0;j<nMaxWFAsignals;j++)
// 	  {
//               cout << "test 1 1 " << j << endl;
// 	      WFA_TimeStructure[0][j] = undefinedValue;
// 	  }
//       }
//       else
//       {
//           cout << "test 1 1a" << endl;
//         cout << "S1 size: " << fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).size() << endl;
// 	  for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).size();j++)
// 	  {
//               cout << "S1 " << j << ": fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(j)" << endl;
// 	      WFA_TimeStructure[0][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS1).at(j);
// 	  }
//       }
//       //S2
//       cout << "test 2" << endl;
//       if (fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS2).size() > nMaxWFAsignals)
//       {
// 	  for (int j=0;j<nMaxWFAsignals;j++)
// 	  {
// 	      WFA_TimeStructure[1][j] = undefinedValue;
// 	  }
//       }
//       else
//       {
// 	  for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS2).size();j++)
// 	  {
// 	      WFA_TimeStructure[1][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS2).at(j);
// 	  }
//       }
//       //S3
//       cout << "test 3" << endl;
//       if (fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS3).size() > nMaxWFAsignals)
//       {
// 	  for (int j=0;j<nMaxWFAsignals;j++)
// 	  {
// 	      WFA_TimeStructure[2][j] = undefinedValue;
// 	  }
//       }
//       else
//       {
// 	  for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS3).size();j++)
// 	  {
// 	      WFA_TimeStructure[2][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eS3).at(j);
// 	  }
//       }
//       //V1
//       cout << "test 4" << endl;
//       if (fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eV1).size() > nMaxWFAsignals)
//       {
// 	  for (int j=0;j<nMaxWFAsignals;j++)
// 	  {
// 	      WFA_TimeStructure[3][j] = undefinedValue;
// 	  }
//       }
//       else
//       {
// 	  for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eV1).size();j++)
// 	  {
// 	      WFA_TimeStructure[3][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eV1).at(j);
// 	  }
//       }
//       //V1p
//       cout << "test 5" << endl;
//       if (fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eV1p).size() > nMaxWFAsignals)
//       {
// 	  for (int j=0;j<nMaxWFAsignals;j++)
// 	  {
// 	      WFA_TimeStructure[4][j] = undefinedValue;
// 	  }
//       }
//       else
//       {
// 	  for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eV1p).size();j++)
// 	  {
// 	      WFA_TimeStructure[4][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::eV1p).at(j);
// 	  }
//       }
//       //PSD
//       cout << "test 6" << endl;
//       if (fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::ePSD).size() > nMaxWFAsignals)
//       {
// 	  for (int j=0;j<nMaxWFAsignals;j++)
// 	  {
// 	      WFA_TimeStructure[5][j] = undefinedValue;
// 	  }
//       }
//       else
//       {
// 	  for (int j=0;j<fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::ePSD).size();j++)
// 	  {
// 	      WFA_TimeStructure[5][j] = fTriggerRaw.GetTimeStructure(det::TimeStructureConst::eWFA,det::TriggerConst::ePSD).at(j);
// 	  }
//       }

      
      
      fTreeQA -> Fill();
  }
  
  void FlowQA::Init_fTreeQA()
  {
    
      //Event
      fTreeQA -> Branch("Event_Id", &Event_Id, "Event_Id/I");
      fTreeQA -> Branch("Run_Id", &Run_Id, "Run_Id/I");
      fTreeQA -> Branch("Event_Timestamp", &Event_Timestamp, "Event_Timestamp/I");
      
//       //SimEvent
// 	//Event
//             fTreeQA -> Branch("Sim_ImpactParameter", &Sim_ImpactParameter, "Sim_ImpactParameter/F");
//             fTreeQA -> Branch("Sim_ReactionPlaneAngle", &Sim_ReactionPlaneAngle, "Sim_ReactionPlaneAngle/F");
//         //Tracks
//             fTreeQA -> Branch("Sim_track_pT", Sim_track_pT, Form("Sim_track_pT[%i]/F",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_phi", Sim_track_phi, Form("Sim_track_phi[%i]/F",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_eta", Sim_track_eta, Form("Sim_track_eta[%i]/F",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_mass", Sim_track_mass, Form("Sim_track_mass[%i]/F",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_charge", Sim_track_charge, Form("Sim_track_charge[%i]/F",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_pdg_id", Sim_track_pdg_id, Form("Sim_track_pdg_id[%i]/I",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_IsInTPC_Total", Sim_track_IsInTPC_Total, Form("Sim_track_IsInTPC_Total[%i]/O",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_IsInTPC_TPCV1", Sim_track_IsInTPC_TPCV1, Form("Sim_track_IsInTPC_TPCV1[%i]/O",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_IsInTPC_TPCV2", Sim_track_IsInTPC_TPCV2, Form("Sim_track_IsInTPC_TPCV2[%i]/O",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_IsInTPC_TPCVmain", Sim_track_IsInTPC_TPCVmain, Form("Sim_track_IsInTPC_TPCVmain[%i]/O",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_Number_of_Hits", Sim_track_Number_of_Hits, Form("Sim_track_Number_of_Hits[%i]/I",nMaxTracks));
//             fTreeQA -> Branch("Sim_track_IsSpectator", Sim_track_IsSpectator, Form("Sim_track_IsSpectator[%i]/O",nMaxTracks));
//             fTreeQA -> Branch("Sim_nTracks", &Sim_nTracks, "Sim_nTracks/I");
//         //PSD
//             fTreeQA -> Branch("Sim_PSD_section_Energy", Sim_PSD_section_Energy, Form("Sim_PSD_section_Energy[%i][%i]/F",nPSD_Modules,nPSD_Sections));
//             fTreeQA -> Branch("Sim_PSD_module_Energy", Sim_PSD_module_Energy, Form("Sim_PSD_module_Energy[%i]/F",nPSD_Modules));
//             fTreeQA -> Branch("Sim_PSD_module_Temperature", Sim_PSD_module_Temperature, Form("Sim_PSD_module_Temperature[%i]/F",nPSD_Modules));
//             fTreeQA -> Branch("Sim_PSD_Energy", &Sim_PSD_Energy, "Sim_PSD_Energy/F");
      
      
      //Main vertex
      fTreeQA -> Branch("Main_Vertex_X", &Main_Vertex_X, "Main_Vertex_X/F");
      fTreeQA -> Branch("Main_Vertex_Y", &Main_Vertex_Y, "Main_Vertex_Y/F");
      fTreeQA -> Branch("Main_Vertex_Z", &Main_Vertex_Z, "Main_Vertex_Z/F");
      
      //Primary vertex
      fTreeQA -> Branch("Primary_Vertex_X", &Primary_Vertex_X, "Primary_Vertex_X/F");
      fTreeQA -> Branch("Primary_Vertex_Y", &Primary_Vertex_Y, "Primary_Vertex_Y/F");
      fTreeQA -> Branch("Primary_Vertex_Z", &Primary_Vertex_Z, "Primary_Vertex_Z/F");
      fTreeQA -> Branch("Has_Primary_Vertex", &Has_Primary_Vertex, "Has_Primary_Vertex/O");
      fTreeQA -> Branch("Primary_Vertex_Q", &Primary_Vertex_Q, "Primary_Vertex_Q/I");
      
      //TPC
      fTreeQA -> Branch("TPC_main_vtx_track_pT", TPC_main_vtx_track_pT, Form("TPC_main_vtx_track_pT[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_px", TPC_main_vtx_track_px, Form("TPC_main_vtx_track_px[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_py", TPC_main_vtx_track_py, Form("TPC_main_vtx_track_py[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_pz", TPC_main_vtx_track_pz, Form("TPC_main_vtx_track_pz[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_eta", TPC_main_vtx_track_eta, Form("TPC_main_vtx_track_eta[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_phi", TPC_main_vtx_track_phi, Form("TPC_main_vtx_track_phi[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_point_x", TPC_main_vtx_track_point_x, Form("TPC_main_vtx_track_point_x[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_point_y", TPC_main_vtx_track_point_y, Form("TPC_main_vtx_track_point_y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_point_z", TPC_main_vtx_track_point_z, Form("TPC_main_vtx_track_point_z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_first_point_x", TPC_main_vtx_track_first_point_x, Form("TPC_main_vtx_track_first_point_x[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_first_point_y", TPC_main_vtx_track_first_point_y, Form("TPC_main_vtx_track_first_point_y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_first_point_z", TPC_main_vtx_track_first_point_z, Form("TPC_main_vtx_track_first_point_z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_last_point_x", TPC_main_vtx_track_last_point_x, Form("TPC_main_vtx_track_last_point_x[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_last_point_y", TPC_main_vtx_track_last_point_y, Form("TPC_main_vtx_track_last_point_y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_last_point_z", TPC_main_vtx_track_last_point_z, Form("TPC_main_vtx_track_last_point_z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_pT", TPC_main_vtx_pT, Form("TPC_main_vtx_pT[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_px", TPC_main_vtx_px, Form("TPC_main_vtx_px[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_py", TPC_main_vtx_py, Form("TPC_main_vtx_py[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_pz", TPC_main_vtx_pz, Form("TPC_main_vtx_pz[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_eta", TPC_main_vtx_eta, Form("TPC_main_vtx_eta[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_phi", TPC_main_vtx_phi, Form("TPC_main_vtx_phi[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_track_px", TPC_main_vtx_ext_track_px, Form("TPC_main_vtx_ext_track_px[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_track_py", TPC_main_vtx_ext_track_py, Form("TPC_main_vtx_ext_track_py[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_track_pz", TPC_main_vtx_ext_track_pz, Form("TPC_main_vtx_ext_track_pz[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_first_track_px", TPC_main_vtx_ext_first_track_px, Form("TPC_main_vtx_ext_first_track_px[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_first_track_py", TPC_main_vtx_ext_first_track_py, Form("TPC_main_vtx_ext_first_track_py[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_first_track_pz", TPC_main_vtx_ext_first_track_pz, Form("TPC_main_vtx_ext_first_track_pz[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_last_track_px", TPC_main_vtx_ext_last_track_px, Form("TPC_main_vtx_ext_last_track_px[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_last_track_py", TPC_main_vtx_ext_last_track_py, Form("TPC_main_vtx_ext_last_track_py[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_last_track_pz", TPC_main_vtx_ext_last_track_pz, Form("TPC_main_vtx_ext_last_track_pz[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_point_x", TPC_main_vtx_ext_point_x, Form("TPC_main_vtx_ext_point_x[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_point_y", TPC_main_vtx_ext_point_y, Form("TPC_main_vtx_ext_point_y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_point_z", TPC_main_vtx_ext_point_z, Form("TPC_main_vtx_ext_point_z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_first_point_x", TPC_main_vtx_ext_first_point_x, Form("TPC_main_vtx_ext_first_point_x[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_first_point_y", TPC_main_vtx_ext_first_point_y, Form("TPC_main_vtx_ext_first_point_y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_first_point_z", TPC_main_vtx_ext_first_point_z, Form("TPC_main_vtx_ext_first_point_z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_last_point_x", TPC_main_vtx_ext_last_point_x, Form("TPC_main_vtx_ext_last_point_x[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_last_point_y", TPC_main_vtx_ext_last_point_y, Form("TPC_main_vtx_ext_last_point_y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_ext_last_point_z", TPC_main_vtx_ext_last_point_z, Form("TPC_main_vtx_ext_last_point_z[%i]/F",nMaxTracks));
	  fTreeQA -> Branch("TPC_main_test_ext_track_px", 		TPC_main_test_ext_track_px,			Form("TPC_main_test_ext_track_px[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_track_py", 		TPC_main_test_ext_track_py,			Form("TPC_main_test_ext_track_py[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_track_pz", 		TPC_main_test_ext_track_pz,			Form("TPC_main_test_ext_track_pz[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_first_track_px", TPC_main_test_ext_first_track_px, 	Form("TPC_main_test_ext_first_track_px[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_first_track_py", TPC_main_test_ext_first_track_py, 	Form("TPC_main_test_ext_first_track_py[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_first_track_pz", TPC_main_test_ext_first_track_pz, 	Form("TPC_main_test_ext_first_track_pz[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_last_track_px", 	TPC_main_test_ext_last_track_px, 	Form("TPC_main_test_ext_last_track_px[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_last_track_py", 	TPC_main_test_ext_last_track_py, 	Form("TPC_main_test_ext_last_track_py[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_last_track_pz", 	TPC_main_test_ext_last_track_pz, 	Form("TPC_main_test_ext_last_track_pz[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_point_x", 		TPC_main_test_ext_point_x, 			Form("TPC_main_test_ext_point_x[%i]/F",			nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_point_y", 		TPC_main_test_ext_point_y, 			Form("TPC_main_test_ext_point_y[%i]/F",			nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_point_z", 		TPC_main_test_ext_point_z, 			Form("TPC_main_test_ext_point_z[%i]/F",			nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_first_point_x", 	TPC_main_test_ext_first_point_x, 	Form("TPC_main_test_ext_first_point_x[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_first_point_y", 	TPC_main_test_ext_first_point_y, 	Form("TPC_main_test_ext_first_point_y[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_first_point_z", 	TPC_main_test_ext_first_point_z, 	Form("TPC_main_test_ext_first_point_z[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_last_point_x", 	TPC_main_test_ext_last_point_x, 	Form("TPC_main_test_ext_last_point_x[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_last_point_y", 	TPC_main_test_ext_last_point_y, 	Form("TPC_main_test_ext_last_point_y[%i]/F",	nMaxTracks));
      fTreeQA -> Branch("TPC_main_test_ext_last_point_z", 	TPC_main_test_ext_last_point_z, 	Form("TPC_main_test_ext_last_point_z[%i]/F",	nMaxTracks));

	  fTreeQA -> Branch("TPC_main_PSD_ext_track_px", 		TPC_main_PSD_ext_track_px,			Form("TPC_main_PSD_ext_track_px[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_PSD_ext_track_py", 		TPC_main_PSD_ext_track_py,			Form("TPC_main_PSD_ext_track_py[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_PSD_ext_track_pz", 		TPC_main_PSD_ext_track_pz,			Form("TPC_main_PSD_ext_track_pz[%i]/F",nMaxTracks));
	  fTreeQA -> Branch("TPC_main_PSD_ext_point_x", 		TPC_main_PSD_ext_point_x, 			Form("TPC_main_PSD_ext_point_x[%i]/F",			nMaxTracks));
      fTreeQA -> Branch("TPC_main_PSD_ext_point_y", 		TPC_main_PSD_ext_point_y, 			Form("TPC_main_PSD_ext_point_y[%i]/F",			nMaxTracks));
      fTreeQA -> Branch("TPC_main_PSD_ext_point_z", 		TPC_main_PSD_ext_point_z, 			Form("TPC_main_PSD_ext_point_z[%i]/F",			nMaxTracks));


      fTreeQA -> Branch("TPC_main_vtx_impact_point_X", TPC_main_vtx_impact_point_X, Form("TPC_main_vtx_impact_point_X[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_impact_point_Y", TPC_main_vtx_impact_point_Y, Form("TPC_main_vtx_impact_point_Y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_impact_point_Z", TPC_main_vtx_impact_point_Z, Form("TPC_main_vtx_impact_point_Z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_pT", TPC_primary_vtx_track_pT, Form("TPC_primary_vtx_track_pT[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_eta", TPC_primary_vtx_track_eta, Form("TPC_primary_vtx_track_eta[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_phi", TPC_primary_vtx_track_phi, Form("TPC_primary_vtx_track_phi[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_pT", TPC_primary_vtx_pT, Form("TPC_primary_vtx_pT[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_eta", TPC_primary_vtx_eta, Form("TPC_primary_vtx_eta[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_phi", TPC_primary_vtx_phi, Form("TPC_primary_vtx_phi[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_impact_point_X", TPC_primary_vtx_impact_point_X, Form("TPC_primary_vtx_impact_point_X[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_impact_point_Y", TPC_primary_vtx_impact_point_Y, Form("TPC_primary_vtx_impact_point_Y[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_impact_point_Z", TPC_primary_vtx_impact_point_Z, Form("TPC_primary_vtx_impact_point_Z[%i]/F",nMaxTracks));
      
      fTreeQA -> Branch("TPC_main_vtx_track_charge", TPC_main_vtx_track_charge, Form("TPC_main_vtx_track_charge[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClusters_Total", TPC_main_vtx_track_nClusters_Total, Form("TPC_main_vtx_track_nClusters_Total[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClusters_TPCV1", TPC_main_vtx_track_nClusters_TPCV1, Form("TPC_main_vtx_track_nClusters_TPCV1[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClusters_TPCV2", TPC_main_vtx_track_nClusters_TPCV2, Form("TPC_main_vtx_track_nClusters_TPCV2[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClusters_TPCVmain", TPC_main_vtx_track_nClusters_TPCVmain, Form("TPC_main_vtx_track_nClusters_TPCVmain[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClusters_TPCVgap", TPC_main_vtx_track_nClusters_TPCVgap, Form("TPC_main_vtx_track_nClusters_TPCVgap[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClustersPotential_Total", TPC_main_vtx_track_nClustersPotential_Total, Form("TPC_main_vtx_track_nClustersPotential_Total[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersPotential_TPCV1", TPC_main_vtx_track_nClustersPotential_TPCV1, Form("TPC_main_vtx_track_nClustersPotential_TPCV1[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersPotential_TPCV2", TPC_main_vtx_track_nClustersPotential_TPCV2, Form("TPC_main_vtx_track_nClustersPotential_TPCV2[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersPotential_TPCVmain", TPC_main_vtx_track_nClustersPotential_TPCVmain, Form("TPC_main_vtx_track_nClustersPotential_TPCVmain[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersPotential_TPCVgap", TPC_main_vtx_track_nClustersPotential_TPCVgap, Form("TPC_main_vtx_track_nClustersPotential_TPCVgap[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersFit_Total", TPC_main_vtx_track_nClustersFit_Total, Form("TPC_main_vtx_track_nClustersFit_Total[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersFit_TPCV1", TPC_main_vtx_track_nClustersFit_TPCV1, Form("TPC_main_vtx_track_nClustersFit_TPCV1[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersFit_TPCV2", TPC_main_vtx_track_nClustersFit_TPCV2, Form("TPC_main_vtx_track_nClustersFit_TPCV2[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersFit_TPCVmain", TPC_main_vtx_track_nClustersFit_TPCVmain, Form("TPC_main_vtx_track_nClustersFit_TPCVmain[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_main_vtx_track_nClustersFit_TPCVgap", TPC_main_vtx_track_nClustersFit_TPCVgap, Form("TPC_main_vtx_track_nClustersFit_TPCVgap[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClustersdEdX_Total", TPC_main_vtx_track_nClustersdEdX_Total, Form("TPC_main_vtx_track_nClustersdEdX_Total[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClustersdEdX_TPCV1", TPC_main_vtx_track_nClustersdEdX_TPCV1, Form("TPC_main_vtx_track_nClustersdEdX_TPCV1[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClustersdEdX_TPCV2", TPC_main_vtx_track_nClustersdEdX_TPCV2, Form("TPC_main_vtx_track_nClustersdEdX_TPCV2[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClustersdEdX_TPCVmain", TPC_main_vtx_track_nClustersdEdX_TPCVmain, Form("TPC_main_vtx_track_nClustersdEdX_TPCVmain[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_nClustersdEdX_TPCVgap", TPC_main_vtx_track_nClustersdEdX_TPCVgap, Form("TPC_main_vtx_track_nClustersdEdX_TPCVgap[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_EnergyClusters_Total", TPC_main_vtx_track_EnergyClusters_Total, Form("TPC_main_vtx_track_EnergyClusters_Total[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_EnergyClusters_TPCV1", TPC_main_vtx_track_EnergyClusters_TPCV1, Form("TPC_main_vtx_track_EnergyClusters_TPCV1[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_EnergyClusters_TPCV2", TPC_main_vtx_track_EnergyClusters_TPCV2, Form("TPC_main_vtx_track_EnergyClusters_TPCV2[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_EnergyClusters_TPCVmain", TPC_main_vtx_track_EnergyClusters_TPCVmain, Form("TPC_main_vtx_track_EnergyClusters_TPCVmain[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_EnergyClusters_TPCVgap", TPC_main_vtx_track_EnergyClusters_TPCVgap, Form("TPC_main_vtx_track_EnergyClusters_TPCVgap[%i]/D",nMaxTracks));      
//       fTreeQA -> Branch("TPC_main_vtx_track_DCAtoVertex_X", TPC_main_vtx_track_DCAtoVertex_X, Form("TPC_main_vtx_track_DCAtoVertex_X[%i]/F",nMaxTracks));   
//       fTreeQA -> Branch("TPC_main_vtx_track_DCAtoVertex_Y", TPC_main_vtx_track_DCAtoVertex_Y, Form("TPC_main_vtx_track_DCAtoVertex_Y[%i]/F",nMaxTracks));   
//       fTreeQA -> Branch("TPC_main_vtx_track_DCAtoVertex_Z", TPC_main_vtx_track_DCAtoVertex_Z, Form("TPC_main_vtx_track_DCAtoVertex_Z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_chi2", TPC_main_vtx_track_chi2, Form("TPC_main_vtx_track_chi2[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_main_vtx_track_ndf", TPC_main_vtx_track_ndf, Form("TPC_main_vtx_track_ndf[%i]/I",nMaxTracks));
      
      fTreeQA -> Branch("TPC_primary_vtx_track_charge", TPC_primary_vtx_track_charge, Form("TPC_primary_vtx_track_charge[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClusters_Total", TPC_primary_vtx_track_nClusters_Total, Form("TPC_primary_vtx_track_nClusters_Total[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClusters_TPCV1", TPC_primary_vtx_track_nClusters_TPCV1, Form("TPC_primary_vtx_track_nClusters_TPCV1[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClusters_TPCV2", TPC_primary_vtx_track_nClusters_TPCV2, Form("TPC_primary_vtx_track_nClusters_TPCV2[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClusters_TPCVmain", TPC_primary_vtx_track_nClusters_TPCVmain, Form("TPC_primary_vtx_track_nClusters_TPCVmain[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClusters_TPCVgap", TPC_primary_vtx_track_nClusters_TPCVgap, Form("TPC_primary_vtx_track_nClusters_TPCVgap[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClustersPotential_Total", TPC_primary_vtx_track_nClustersPotential_Total, Form("TPC_primary_vtx_track_nClustersPotential_Total[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersPotential_TPCV1", TPC_primary_vtx_track_nClustersPotential_TPCV1, Form("TPC_primary_vtx_track_nClustersPotential_TPCV1[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersPotential_TPCV2", TPC_primary_vtx_track_nClustersPotential_TPCV2, Form("TPC_primary_vtx_track_nClustersPotential_TPCV2[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersPotential_TPCVmain", TPC_primary_vtx_track_nClustersPotential_TPCVmain, Form("TPC_primary_vtx_track_nClustersPotential_TPCVmain[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersPotential_TPCVgap", TPC_primary_vtx_track_nClustersPotential_TPCVgap, Form("TPC_primary_vtx_track_nClustersPotential_TPCVgap[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersFit_Total", TPC_primary_vtx_track_nClustersFit_Total, Form("TPC_primary_vtx_track_nClustersFit_Total[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersFit_TPCV1", TPC_primary_vtx_track_nClustersFit_TPCV1, Form("TPC_primary_vtx_track_nClustersFit_TPCV1[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersFit_TPCV2", TPC_primary_vtx_track_nClustersFit_TPCV2, Form("TPC_primary_vtx_track_nClustersFit_TPCV2[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersFit_TPCVmain", TPC_primary_vtx_track_nClustersFit_TPCVmain, Form("TPC_primary_vtx_track_nClustersFit_TPCVmain[%i]/I",nMaxTracks));
//       fTreeQA -> Branch("TPC_primary_vtx_track_nClustersFit_TPCVgap", TPC_primary_vtx_track_nClustersFit_TPCVgap, Form("TPC_primary_vtx_track_nClustersFit_TPCVgap[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClustersdEdX_Total", TPC_primary_vtx_track_nClustersdEdX_Total, Form("TPC_primary_vtx_track_nClustersdEdX_Total[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClustersdEdX_TPCV1", TPC_primary_vtx_track_nClustersdEdX_TPCV1, Form("TPC_primary_vtx_track_nClustersdEdX_TPCV1[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClustersdEdX_TPCV2", TPC_primary_vtx_track_nClustersdEdX_TPCV2, Form("TPC_primary_vtx_track_nClustersdEdX_TPCV2[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClustersdEdX_TPCVmain", TPC_primary_vtx_track_nClustersdEdX_TPCVmain, Form("TPC_primary_vtx_track_nClustersdEdX_TPCVmain[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_nClustersdEdX_TPCVgap", TPC_primary_vtx_track_nClustersdEdX_TPCVgap, Form("TPC_primary_vtx_track_nClustersdEdX_TPCVgap[%i]/I",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_EnergyClusters_Total", TPC_primary_vtx_track_EnergyClusters_Total, Form("TPC_primary_vtx_track_EnergyClusters_Total[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_EnergyClusters_TPCV1", TPC_primary_vtx_track_EnergyClusters_TPCV1, Form("TPC_primary_vtx_track_EnergyClusters_TPCV1[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_EnergyClusters_TPCV2", TPC_primary_vtx_track_EnergyClusters_TPCV2, Form("TPC_primary_vtx_track_EnergyClusters_TPCV2[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_EnergyClusters_TPCVmain", TPC_primary_vtx_track_EnergyClusters_TPCVmain, Form("TPC_primary_vtx_track_EnergyClusters_TPCVmain[%i]/D",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_EnergyClusters_TPCVgap", TPC_primary_vtx_track_EnergyClusters_TPCVgap, Form("TPC_primary_vtx_track_EnergyClusters_TPCVgap[%i]/D",nMaxTracks));      
//       fTreeQA -> Branch("TPC_primary_vtx_track_DCAtoVertex_X", TPC_primary_vtx_track_DCAtoVertex_X, Form("TPC_primary_vtx_track_DCAtoVertex_X[%i]/F",nMaxTracks));   
//       fTreeQA -> Branch("TPC_primary_vtx_track_DCAtoVertex_Y", TPC_primary_vtx_track_DCAtoVertex_Y, Form("TPC_primary_vtx_track_DCAtoVertex_Y[%i]/F",nMaxTracks));   
//       fTreeQA -> Branch("TPC_primary_vtx_track_DCAtoVertex_Z", TPC_primary_vtx_track_DCAtoVertex_Z, Form("TPC_primary_vtx_track_DCAtoVertex_Z[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_chi2", TPC_primary_vtx_track_chi2, Form("TPC_primary_vtx_track_chi2[%i]/F",nMaxTracks));
      fTreeQA -> Branch("TPC_primary_vtx_track_ndf", TPC_primary_vtx_track_ndf, Form("TPC_primary_vtx_track_ndf[%i]/I",nMaxTracks));
      
      
      fTreeQA -> Branch("TPC_main_vtx_nTracks", &TPC_main_vtx_nTracks, "TPC_main_vtx_nTracks/I");
      fTreeQA -> Branch("TPC_primary_vtx_nTracks", &TPC_primary_vtx_nTracks, "TPC_primary_vtx_nTracks/I");
//       fTreeQA -> Branch("TPC_Multiplicity", &TPC_Multiplicity, "TPC_Multiplicity/I");
//       fTreeQA -> Branch("TPC_Multiplicity_all", &TPC_Multiplicity_all, "TPC_Multiplicity_all/I");
//       fTreeQA -> Branch("Primary_Multiplicity_all", &Primary_Multiplicity_all, "Primary_Multiplicity_all/I");
//       fTreeQA -> Branch("TPC_Multiplicity_Clusters_VTPC1_VTPC2", &TPC_Multiplicity_Clusters_VTPC1_VTPC2, "TPC_Multiplicity_Clusters_VTPC1_VTPC2/I");
//       fTreeQA -> Branch("TPC_Multiplicity_Clusters_All", &TPC_Multiplicity_Clusters_All, "TPC_Multiplicity_Clusters_All/I");
      
//       fTreeQA -> Branch("TPC_cos1", &TPC_cos1, "TPC_cos1/F");
//       fTreeQA -> Branch("TPC_sin1", &TPC_sin1, "TPC_sin1/F");
//       fTreeQA -> Branch("TPC_cos2", &TPC_cos2, "TPC_cos2/F");
//       fTreeQA -> Branch("TPC_sin2", &TPC_sin2, "TPC_sin2/F");
      
      
      //PSD
      fTreeQA -> Branch("PSD_module_Number", &PSD_module_Number, "PSD_module_Number/I");
      fTreeQA -> Branch("PSD_section_Number", &PSD_section_Number,"PSD_section_Number/I");
//       fTreeQA -> Branch("PSD_section_slice_Energy", PSD_section_slice_Energy, Form("PSD_section_slice_Energy[%i]/F",nPSD_Sections));
      fTreeQA -> Branch("PSD_module_X", PSD_module_X, Form("PSD_module_X[%i]/F",nPSD_Modules));
      fTreeQA -> Branch("PSD_module_Y", PSD_module_Y, Form("PSD_module_Y[%i]/F",nPSD_Modules));
      fTreeQA -> Branch("PSD_module_Z", PSD_module_Z, Form("PSD_module_Z[%i]/F",nPSD_Modules));
      fTreeQA -> Branch("PSD_module_Energy", PSD_module_Energy, Form("PSD_module_Energy[%i]/F",nPSD_Modules));
      fTreeQA -> Branch("PSD_module_Energy_default", PSD_module_Energy_default, Form("PSD_module_Energy_default[%i]/F",nPSD_Modules));
      fTreeQA -> Branch("PSD_module_number_of_sections", PSD_module_number_of_sections, Form("PSD_module_number_of_sections[%i]/I",nPSD_Modules));
      fTreeQA -> Branch("PSD_section_Energy", PSD_section_Energy, Form("PSD_section_Energy[%i][%i]/F",nPSD_Modules,nPSD_Sections));
      fTreeQA -> Branch("PSD_section_Number_array", PSD_section_Number_array, Form("PSD_section_Number_array[%i][%i]/I",nPSD_Modules,nPSD_Sections));
      fTreeQA -> Branch("PSD_Energy", &PSD_Energy, "PSD_Energy/F");
      fTreeQA -> Branch("PSD_1_Energy", &PSD_1_Energy, "PSD_1_Energy/F");
      fTreeQA -> Branch("PSD_2_Energy", &PSD_2_Energy, "PSD_2_Energy/F");
      fTreeQA -> Branch("PSD_3_Energy", &PSD_3_Energy, "PSD_3_Energy/F");
      
      //Primary vertex
//       fTreeQA -> Branch("Vertex_X", &Vertex_X, "Vertex_X/F");
//       fTreeQA -> Branch("Vertex_Y", &Vertex_Y, "Vertex_Y/F");
//       fTreeQA -> Branch("Vertex_Z", &Vertex_Z, "Vertex_Z/F");
      
//       //Triggers (old style variables)
//       fTreeQA -> Branch("T1", &T1, "T1/O");
//       fTreeQA -> Branch("T2", &T2, "T2/O");
//       fTreeQA -> Branch("T3", &T3, "T3/O");
//       fTreeQA -> Branch("T4", &T4, "T4/O");
      
      //Triggers
      fTreeQA -> Branch("BPD_Status", BPD_Status, Form("BPD_Status[%i][%i][%i]/I",nRawReco,nBPD,nBPDcomponents));  
      fTreeQA -> Branch("BPD_Position", BPD_Position, Form("BPD_Position[%i][%i][%i]/F",nRawReco,nBPD,nBPDcomponents)); 
      fTreeQA -> Branch("BPD_PositionError", BPD_PositionError, Form("BPD_PositionError[%i][%i][%i]/F",nRawReco,nBPD,nBPDcomponents)); 
//       fTreeQA -> Branch("BPD_Z", BPD_Z, Form("BPD_Z[%i][%i][%i]/F",nRawReco,nBPD,nBPDcomponents)); 
//       fTreeQA -> Branch("BPD_RMS", BPD_RMS, Form("BPD_RMS[%i][%i][%i]/F",nRawReco,nBPD,nBPDcomponents)); 
//       fTreeQA -> Branch("BPD_Maximum", BPD_Maximum, Form("BPD_Maximum[%i][%i][%i]/F",nRawReco,nBPD,nBPDcomponents)); 
//       fTreeQA -> Branch("BPD_Charge", BPD_Charge, Form("BPD_Charge[%i][%i][%i]/F",nRawReco,nBPD,nBPDcomponents)); 
//       fTreeQA -> Branch("BPD_SumOfAll", BPD_SumOfAll, Form("BPD_SumOfAll[%i][%i][%i]/F",nRawReco,nBPD,nBPDcomponents)); 
      fTreeQA -> Branch("triggersADC", triggersADC, Form("triggersADC[%i][%i]/F",nRawReco,nTriggers_Simple)); 
      fTreeQA -> Branch("isTriggers_Simple", isTriggers_Simple, Form("isTriggers_Simple[%i][%i]/O",nRawReco,nTriggers_Simple)); 
      fTreeQA -> Branch("isTriggers_Combined", isTriggers_Combined, Form("isTriggers_Combined[%i][%i]/O",nRawReco,nTriggers_Combined)); 
//       fTreeQA -> Branch("Beam_Momentum", Beam_Momentum, Form("Beam_Momentum[%i][%i]/F",nRawReco,nBeamComponents)); 
//       fTreeQA -> Branch("Beam_Fitted2DLineXZ", Beam_Fitted2DLineXZ, Form("Beam_Fitted2DLineXZ[%i][%i]/F",nRawReco,nBeamComponents)); 
//       fTreeQA -> Branch("Beam_Fitted2DLineYZ", Beam_Fitted2DLineYZ, Form("Beam_Fitted2DLineYZ[%i][%i]/F",nRawReco,nBeamComponents)); 
//       fTreeQA -> Branch("Beam_Status", Beam_Status, Form("Beam_Status[%i]/I",nRawReco)); 
      fTreeQA -> Branch("WFA_TimeStructure", WFA_TimeStructure, Form("WFA_TimeStructure[%i][%i]/F",nTriggers_Simple,nMaxWFAsignals)); 
      fTreeQA -> Branch("WFA_NumberOfSignalHits", WFA_NumberOfSignalHits, Form("WFA_NumberOfSignalHits[%i]/I",nTriggers_Simple)); 
      fTreeQA -> Branch("FitVertexX", &FitVertexX, "FitVertexX/F");
      fTreeQA -> Branch("FitVertexY", &FitVertexY, "FitVertexY/F");
      fTreeQA -> Branch("FitVertexZ", &FitVertexZ, "FitVertexZ/F");
      fTreeQA -> Branch("FitVertexQ", &FitVertexQ, "FitVertexQ/I");
  }
  
  void FlowQA::ClearEvent()
  {
//       fTreeQA->ResetBranchAddresses();
      Event_Timestamp = 0;
      Event_Id = 0;
      Run_Id = 0;
      
      for (int i=0;i<nMaxTracks;i++)
      {
          
           Sim_track_pT[i] = undefinedValue;
           Sim_track_phi[i] = undefinedValue;
           Sim_track_eta[i] = undefinedValue;
           Sim_track_mass[i] = undefinedValue;
           Sim_track_charge[i] = undefinedValue;
           Sim_track_pdg_id[i] = undefinedValue;
           Sim_track_IsInTPC_Total[i] = undefinedValue;
           Sim_track_IsInTPC_TPCV1[i] = undefinedValue;
           Sim_track_IsInTPC_TPCV2[i] = undefinedValue;
           Sim_track_IsInTPC_TPCVmain[i] = undefinedValue;
           Sim_track_Number_of_Hits[i] = undefinedValue;
           Sim_track_IsSpectator[i] = undefinedValue;
	    
	  TPC_main_vtx_track_pT[i] = undefinedValue;
      TPC_main_vtx_track_px[i] = undefinedValue;
      TPC_main_vtx_track_py[i] = undefinedValue;
      TPC_main_vtx_track_pz[i] = undefinedValue;
	  TPC_main_vtx_track_eta[i] = undefinedValue;
	  TPC_main_vtx_track_phi[i] = undefinedValue;
      TPC_main_vtx_track_point_x[i] = undefinedValue;
      TPC_main_vtx_track_point_y[i] = undefinedValue;
      TPC_main_vtx_track_point_z[i] = undefinedValue;
      TPC_main_vtx_track_first_point_x[i] = undefinedValue;
      TPC_main_vtx_track_first_point_y[i] = undefinedValue;
      TPC_main_vtx_track_first_point_z[i] = undefinedValue;
      TPC_main_vtx_track_last_point_x[i] = undefinedValue;
      TPC_main_vtx_track_last_point_y[i] = undefinedValue;
      TPC_main_vtx_track_last_point_z[i] = undefinedValue;

	  TPC_main_vtx_pT[i] = undefinedValue;
      TPC_main_vtx_px[i] = undefinedValue;
      TPC_main_vtx_py[i] = undefinedValue;
      TPC_main_vtx_pz[i] = undefinedValue;
	  TPC_main_vtx_eta[i] = undefinedValue;
	  TPC_main_vtx_phi[i] = undefinedValue;
	  TPC_main_vtx_impact_point_X[i] = undefinedValue;
	  TPC_main_vtx_impact_point_Y[i] = undefinedValue;
	  TPC_main_vtx_impact_point_Z[i] = undefinedValue;

      TPC_main_vtx_ext_track_px[i] = undefinedValue;
      TPC_main_vtx_ext_track_py[i] = undefinedValue;
      TPC_main_vtx_ext_track_pz[i] = undefinedValue;
      TPC_main_vtx_ext_first_track_px[i] = undefinedValue;
      TPC_main_vtx_ext_first_track_py[i] = undefinedValue;
      TPC_main_vtx_ext_first_track_pz[i] = undefinedValue;
      TPC_main_vtx_ext_last_track_px[i] = undefinedValue;
      TPC_main_vtx_ext_last_track_py[i] = undefinedValue;
      TPC_main_vtx_ext_last_track_pz[i] = undefinedValue;
      TPC_main_vtx_ext_point_x[i] = undefinedValue;
      TPC_main_vtx_ext_point_y[i] = undefinedValue;
      TPC_main_vtx_ext_point_z[i] = undefinedValue;
      TPC_main_vtx_ext_first_point_x[i] = undefinedValue;
      TPC_main_vtx_ext_first_point_y[i] = undefinedValue;
      TPC_main_vtx_ext_first_point_z[i] = undefinedValue;
      TPC_main_vtx_ext_last_point_x[i] = undefinedValue;
      TPC_main_vtx_ext_last_point_y[i] = undefinedValue;
      TPC_main_vtx_ext_last_point_z[i] = undefinedValue;

	  TPC_main_test_ext_track_px[i] = undefinedValue;
      TPC_main_test_ext_track_py[i] = undefinedValue;
      TPC_main_test_ext_track_pz[i] = undefinedValue;
      TPC_main_test_ext_first_track_px[i] = undefinedValue;
      TPC_main_test_ext_first_track_py[i] = undefinedValue;
      TPC_main_test_ext_first_track_pz[i] = undefinedValue;
      TPC_main_test_ext_last_track_px[i] = undefinedValue;
      TPC_main_test_ext_last_track_py[i] = undefinedValue;
      TPC_main_test_ext_last_track_pz[i] = undefinedValue;
      TPC_main_test_ext_point_x[i] = undefinedValue;
      TPC_main_test_ext_point_y[i] = undefinedValue;
      TPC_main_test_ext_point_z[i] = undefinedValue;
      TPC_main_test_ext_first_point_x[i] = undefinedValue;
      TPC_main_test_ext_first_point_y[i] = undefinedValue;
      TPC_main_test_ext_first_point_z[i] = undefinedValue;
      TPC_main_test_ext_last_point_x[i] = undefinedValue;
      TPC_main_test_ext_last_point_y[i] = undefinedValue;
      TPC_main_test_ext_last_point_z[i] = undefinedValue;

	  TPC_main_PSD_ext_track_px[i] = undefinedValue;
      TPC_main_PSD_ext_track_py[i] = undefinedValue;
      TPC_main_PSD_ext_track_pz[i] = undefinedValue;
	  TPC_main_PSD_ext_point_x[i] = undefinedValue;
      TPC_main_PSD_ext_point_y[i] = undefinedValue;
      TPC_main_PSD_ext_point_z[i] = undefinedValue;


	  TPC_primary_vtx_track_pT[i] = undefinedValue;
	  TPC_primary_vtx_track_eta[i] = undefinedValue;
	  TPC_primary_vtx_track_phi[i] = undefinedValue;

	  TPC_primary_vtx_pT[i] = undefinedValue;
	  TPC_primary_vtx_eta[i] = undefinedValue;
	  TPC_primary_vtx_phi[i] = undefinedValue;
	  TPC_primary_vtx_impact_point_X[i] = undefinedValue;
	  TPC_primary_vtx_impact_point_Y[i] = undefinedValue;
	  TPC_primary_vtx_impact_point_Z[i] = undefinedValue;
	
	   TPC_main_vtx_track_charge[i] = undefinedValue;
	   TPC_main_vtx_track_nClusters_Total[i] = undefinedValue;
	   TPC_main_vtx_track_nClusters_TPCV1[i] = undefinedValue;
	   TPC_main_vtx_track_nClusters_TPCV2[i] = undefinedValue;
	   TPC_main_vtx_track_nClusters_TPCVmain[i] = undefinedValue;
	   TPC_main_vtx_track_nClusters_TPCVgap[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersPotential_Total[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersPotential_TPCV1[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersPotential_TPCV2[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersPotential_TPCVmain[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersPotential_TPCVgap[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersFit_Total[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersFit_TPCV1[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersFit_TPCV2[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersFit_TPCVmain[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersFit_TPCVgap[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersdEdX_Total[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersdEdX_TPCV1[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersdEdX_TPCV2[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersdEdX_TPCVmain[i] = undefinedValue;
	   TPC_main_vtx_track_nClustersdEdX_TPCVgap[i] = undefinedValue;
	   TPC_main_vtx_track_EnergyClusters_Total[i] = undefinedValue;
	   TPC_main_vtx_track_EnergyClusters_TPCV1[i] = undefinedValue;
	   TPC_main_vtx_track_EnergyClusters_TPCV2[i] = undefinedValue;
	   TPC_main_vtx_track_EnergyClusters_TPCVmain[i] = undefinedValue;
	   TPC_main_vtx_track_EnergyClusters_TPCVgap[i] = undefinedValue;
	   TPC_main_vtx_track_DCAtoVertex_X[i] = undefinedValue;
	   TPC_main_vtx_track_DCAtoVertex_Y[i] = undefinedValue;
	   TPC_main_vtx_track_DCAtoVertex_Z[i] = undefinedValue;
	   TPC_main_vtx_track_chi2[i] = undefinedValue;  
	   TPC_main_vtx_track_ndf[i] = undefinedValue;  
	   
	   TPC_primary_vtx_track_charge[i] = undefinedValue;
	   TPC_primary_vtx_track_nClusters_Total[i] = undefinedValue;
	   TPC_primary_vtx_track_nClusters_TPCV1[i] = undefinedValue;
	   TPC_primary_vtx_track_nClusters_TPCV2[i] = undefinedValue;
	   TPC_primary_vtx_track_nClusters_TPCVmain[i] = undefinedValue;
	   TPC_primary_vtx_track_nClusters_TPCVgap[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersPotential_Total[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersPotential_TPCV1[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersPotential_TPCV2[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersPotential_TPCVmain[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersPotential_TPCVgap[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersFit_Total[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersFit_TPCV1[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersFit_TPCV2[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersFit_TPCVmain[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersFit_TPCVgap[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersdEdX_Total[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersdEdX_TPCV1[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersdEdX_TPCV2[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersdEdX_TPCVmain[i] = undefinedValue;
	   TPC_primary_vtx_track_nClustersdEdX_TPCVgap[i] = undefinedValue;
	   TPC_primary_vtx_track_EnergyClusters_Total[i] = undefinedValue;
	   TPC_primary_vtx_track_EnergyClusters_TPCV1[i] = undefinedValue;
	   TPC_primary_vtx_track_EnergyClusters_TPCV2[i] = undefinedValue;
	   TPC_primary_vtx_track_EnergyClusters_TPCVmain[i] = undefinedValue;
	   TPC_primary_vtx_track_EnergyClusters_TPCVgap[i] = undefinedValue;
	   TPC_primary_vtx_track_DCAtoVertex_X[i] = undefinedValue;
	   TPC_primary_vtx_track_DCAtoVertex_Y[i] = undefinedValue;
	   TPC_primary_vtx_track_DCAtoVertex_Z[i] = undefinedValue;
	   TPC_primary_vtx_track_chi2[i] = undefinedValue;  
	   TPC_primary_vtx_track_ndf[i] = undefinedValue;  
      }
      Sim_nTracks = 0;
      Sim_PSD_Energy = 0;
      TPC_main_vtx_nTracks = 0;
      TPC_primary_vtx_nTracks = 0;
      TPC_Multiplicity = 0;
      TPC_Multiplicity_all = 0;
      Primary_Multiplicity_all = 0;
      TPC_Multiplicity_Clusters_All = 0;
      TPC_Multiplicity_Clusters_VTPC1_VTPC2 = 0;
      TPC_cos1=0;
      TPC_sin1=0;
      TPC_cos2=0;
      TPC_sin2=0;
      
      PSD_Energy = 0;
      PSD_1_Energy = 0;
      PSD_2_Energy = 0;
      PSD_3_Energy = 0;
      PSD_section_Number = 0;
      PSD_module_Number = 0;
      for (int i=0;i<nPSD_Modules;i++)
      {
	  for (int j=0;j<nPSD_Sections;j++)
	  {
	      PSD_section_Energy[i][j] = 0;
	      PSD_section_Number_array[i][j] = 0;
              Sim_PSD_section_Energy[i][j] = 0;
	  }
	  PSD_module_X[i] = 0;
	  PSD_module_Y[i] = 0;
	  PSD_module_Z[i] = 17000;
          PSD_module_Energy[i] = 0;
	  PSD_module_Energy_default[i] = 0;
	  PSD_module_number_of_sections[i] = 0;
          Sim_PSD_module_Energy[i] = 0;
          Sim_PSD_module_Temperature[i] = 0;
      }
      PSD_module_Z[44] = 17000.0-259.8;
      for (int i=0;i<nPSD_Sections;i++)
      {
	  PSD_section_slice_Energy[i] = 0;
      }
      
      T1 = false;
      T2 = false;
      T3 = false;
      T4 = false;
      
      for (int k=0;k<nRawReco;k++)
      {
	  for (int i=0;i<nBPD;i++)
	  {
	      for (int j=0;j<nBPDcomponents;j++)
	      {
		  BPD_Status[k][i][j] = undefinedValue;
		  BPD_Position[k][i][j] = undefinedValue;
		  BPD_PositionError[k][i][j] = undefinedValue;
		  BPD_Z[k][i][j] = undefinedValue;
		  BPD_RMS[k][i][j] = undefinedValue;
		  BPD_Maximum[k][i][j] = undefinedValue;
		  BPD_Charge[k][i][j] = undefinedValue;
		  BPD_SumOfAll[k][i][j] = undefinedValue;
	      }
	  }
	  for (int i=0;i<nTriggers_Simple;i++)
	  {
	      triggersADC[k][i] = undefinedValue;
	      isTriggers_Simple[k][i] = false;
	  }
	  for (int i=0;i<nTriggers_Combined;i++)
	  {
	      isTriggers_Combined[k][i] = false;
	  }
	  for (int i=0;i<nBeamComponents;i++)
	  {
	      Beam_Momentum[k][i] = undefinedValue;
	      Beam_Fitted2DLineXZ[k][i] = undefinedValue;
	      Beam_Fitted2DLineYZ[k][i] = undefinedValue;
	  }
	  Beam_Status[k] = undefinedValue;
      }
      for (int i=0;i<nTriggers_Simple;i++)
      {
	  for (int j=0;j<nMaxWFAsignals;j++)
	  {
	      WFA_TimeStructure[i][j] = undefinedValue;
	  }
	  WFA_NumberOfSignalHits[i] = undefinedValue;
      }
      FitVertexX = undefinedValue;
      FitVertexY = undefinedValue;
      FitVertexZ = undefinedValue;
      FitVertexQ = undefinedValue;  
  }
  
  
}
