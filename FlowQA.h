/**

   declaration for FlowQA
   date 22/10/2012

*/

#ifndef _FlowQA_FlowQA_h_
#define _FlowQA_FlowQA_h_

#include <fwk/VModule.h>
#include <utl/Branch.h>
#include <map>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <evt/rec/Cluster.h>
#include <TNtuple.h>
#include <TMath.h> //vblinov
#include <TTree.h> //vblinov
#include "bitset"
#include <vector>

namespace FlowQA {

  /**
     \class FlowQA
     \author V. Blinov
     \brief Data quality check for Flow Performance studies
     \ingroup MonitoringModules
  */

  class FlowQA : public fwk::VModule {

  public:
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();

  private:

    enum QA_Tree {
      nMaxTracks = 2000,
      undefinedValue = -999,
      nPSD_Sections = 10,
      nPSD_Modules = 45,
      nBPDcomponents = 3, 		//0=x, 1=y, 2=z,
      nTriggers_Simple = 6, 		// s1, s2, s3, v1, v1p, PSD
      nTriggers_Combined = 4, 		// T1, T2, T3, T4
      nBPD = 3, 			//BPD1, BPD2, BPD3
      nRawReco = 2, 			//0=raw, 1=reco
      nBeamComponents = 3, 		//0=x, 1=y, 2=z
      nMaxWFAsignals = 2000
    };
    enum trigs{
		eTrig1,
		eTrig2,
		eTrig3,
		eTrig4
    };
       
  private:
    void Calculate_QA(evt::Event& event);
    void Init_fTreeQA();
    void ClearEvent();
    TTree* fTreeQA;
    
    //init
    std::string fQAFileName;
    double eBeam;
    int layerNum;
    int maxpad;
    //Counters
    std::bitset<8> trigBit;
    Bool_t timeUsed;
    unsigned int timeStart;
    int timeStamp;
    unsigned int fNEvents;

    unsigned int fNT1;
    unsigned int fNT2;
    unsigned int fNT3;
    unsigned int plotDriftVel;


    unsigned int fNGoodBPD;
    unsigned int fNHasMainVertex;
    unsigned int fNHasFittedVertex;
    unsigned int fNMainVertexPosition;
    unsigned int fNVertexTracks;

    unsigned int fNTracks;
    unsigned int fNTrackStatus;
    unsigned int fNVTPCClusters;
    unsigned int fNTotClusters;
    unsigned int fNClusterRatio;
    unsigned int fNImpactParameter;
    unsigned int fNHasTOF;
    
    // --- cut values
    bool fEventCutsOn;
    bool fUseVertexTracks;
    
    bool fLoopThroughMainVertex;
    unsigned int fNMinVertexTracks;
    unsigned int fNMinVTPCClusters;
    unsigned int fNMinTotalClusters;
    double fMinPointRatio;
    double fMaxImpactX;
    double fMaxImpactY;
    double fZMin;
    double fZMax;
    
    double dVel[5];
    
    
    //Event
    int Event_Timestamp;
    int Event_Id;
    int Run_Id;
    
    //SimEvent
	//Event
	float Sim_ImpactParameter;
	float Sim_ReactionPlaneAngle;
        //Tracks
        float Sim_track_pT[nMaxTracks];
        float Sim_track_phi[nMaxTracks];
        float Sim_track_eta[nMaxTracks];
        float Sim_track_mass[nMaxTracks];
        float Sim_track_charge[nMaxTracks];
        int Sim_track_pdg_id[nMaxTracks];
        bool Sim_track_IsInTPC_Total[nMaxTracks];
        bool Sim_track_IsInTPC_TPCV1[nMaxTracks];
        bool Sim_track_IsInTPC_TPCV2[nMaxTracks];
        bool Sim_track_IsInTPC_TPCVmain[nMaxTracks];
        int Sim_track_Number_of_Hits[nMaxTracks];
        bool Sim_track_IsSpectator[nMaxTracks];
        int Sim_nTracks;
        
        //PSD
        float Sim_PSD_section_Energy[nPSD_Modules][nPSD_Sections];
        float Sim_PSD_module_Energy[nPSD_Modules];
        float Sim_PSD_module_Temperature[nPSD_Modules];
        float Sim_PSD_Energy;
    //RecEvent    
        //TPC
        int TPC_main_vtx_nTracks;
        int TPC_primary_vtx_nTracks;
        
        float TPC_main_vtx_track_pT[nMaxTracks];
        float TPC_main_vtx_track_px[nMaxTracks];
        float TPC_main_vtx_track_py[nMaxTracks];
        float TPC_main_vtx_track_pz[nMaxTracks];
        float TPC_main_vtx_track_eta[nMaxTracks];
        float TPC_main_vtx_track_phi[nMaxTracks];
        float TPC_main_vtx_track_point_x[nMaxTracks];
        float TPC_main_vtx_track_point_y[nMaxTracks];
        float TPC_main_vtx_track_point_z[nMaxTracks];
        float TPC_main_vtx_track_first_point_x[nMaxTracks];
        float TPC_main_vtx_track_first_point_y[nMaxTracks];
        float TPC_main_vtx_track_first_point_z[nMaxTracks];
        float TPC_main_vtx_track_last_point_x[nMaxTracks];
        float TPC_main_vtx_track_last_point_y[nMaxTracks];
        float TPC_main_vtx_track_last_point_z[nMaxTracks];
        float TPC_main_vtx_ext_track_px[nMaxTracks];
        float TPC_main_vtx_ext_track_py[nMaxTracks];
        float TPC_main_vtx_ext_track_pz[nMaxTracks];
        float TPC_main_vtx_ext_first_track_px[nMaxTracks];
        float TPC_main_vtx_ext_first_track_py[nMaxTracks];
        float TPC_main_vtx_ext_first_track_pz[nMaxTracks];
        float TPC_main_vtx_ext_last_track_px[nMaxTracks];
        float TPC_main_vtx_ext_last_track_py[nMaxTracks];
        float TPC_main_vtx_ext_last_track_pz[nMaxTracks];
        float TPC_main_vtx_ext_point_x[nMaxTracks];
        float TPC_main_vtx_ext_point_y[nMaxTracks];
        float TPC_main_vtx_ext_point_z[nMaxTracks];
        float TPC_main_vtx_ext_first_point_x[nMaxTracks];
        float TPC_main_vtx_ext_first_point_y[nMaxTracks];
        float TPC_main_vtx_ext_first_point_z[nMaxTracks];
        float TPC_main_vtx_ext_last_point_x[nMaxTracks];
        float TPC_main_vtx_ext_last_point_y[nMaxTracks];
        float TPC_main_vtx_ext_last_point_z[nMaxTracks];

		float TPC_main_test_ext_track_px[nMaxTracks];
        float TPC_main_test_ext_track_py[nMaxTracks];
        float TPC_main_test_ext_track_pz[nMaxTracks];
        float TPC_main_test_ext_first_track_px[nMaxTracks];
        float TPC_main_test_ext_first_track_py[nMaxTracks];
        float TPC_main_test_ext_first_track_pz[nMaxTracks];
        float TPC_main_test_ext_last_track_px[nMaxTracks];
        float TPC_main_test_ext_last_track_py[nMaxTracks];
        float TPC_main_test_ext_last_track_pz[nMaxTracks];
        float TPC_main_test_ext_point_x[nMaxTracks];
        float TPC_main_test_ext_point_y[nMaxTracks];
        float TPC_main_test_ext_point_z[nMaxTracks];
        float TPC_main_test_ext_first_point_x[nMaxTracks];
        float TPC_main_test_ext_first_point_y[nMaxTracks];
        float TPC_main_test_ext_first_point_z[nMaxTracks];
        float TPC_main_test_ext_last_point_x[nMaxTracks];
        float TPC_main_test_ext_last_point_y[nMaxTracks];
        float TPC_main_test_ext_last_point_z[nMaxTracks];

		float TPC_main_PSD_ext_track_px[nMaxTracks];
        float TPC_main_PSD_ext_track_py[nMaxTracks];
        float TPC_main_PSD_ext_track_pz[nMaxTracks];
		float TPC_main_PSD_ext_point_x[nMaxTracks];
        float TPC_main_PSD_ext_point_y[nMaxTracks];
        float TPC_main_PSD_ext_point_z[nMaxTracks];
	
        float TPC_main_vtx_pT[nMaxTracks];
        float TPC_main_vtx_px[nMaxTracks];
        float TPC_main_vtx_py[nMaxTracks];
        float TPC_main_vtx_pz[nMaxTracks];
        float TPC_main_vtx_eta[nMaxTracks];
        float TPC_main_vtx_phi[nMaxTracks];
	float TPC_main_vtx_impact_point_X[nMaxTracks];
	float TPC_main_vtx_impact_point_Y[nMaxTracks];
	float TPC_main_vtx_impact_point_Z[nMaxTracks];
	
        float TPC_primary_vtx_track_pT[nMaxTracks];
        float TPC_primary_vtx_track_eta[nMaxTracks];
        float TPC_primary_vtx_track_phi[nMaxTracks];
	
        float TPC_primary_vtx_pT[nMaxTracks];
        float TPC_primary_vtx_eta[nMaxTracks];
        float TPC_primary_vtx_phi[nMaxTracks];
	float TPC_primary_vtx_impact_point_X[nMaxTracks];
	float TPC_primary_vtx_impact_point_Y[nMaxTracks];
	float TPC_primary_vtx_impact_point_Z[nMaxTracks];
	
        int TPC_main_vtx_track_charge[nMaxTracks];
        double TPC_main_vtx_track_dEdx[nMaxTracks];
        int TPC_main_vtx_track_nClusters_Total[nMaxTracks];
        int TPC_main_vtx_track_nClusters_TPCV1[nMaxTracks];
        int TPC_main_vtx_track_nClusters_TPCV2[nMaxTracks];
        int TPC_main_vtx_track_nClusters_TPCVmain[nMaxTracks];
        int TPC_main_vtx_track_nClusters_TPCVgap[nMaxTracks];
        int TPC_main_vtx_track_nClustersPotential_Total[nMaxTracks];
        int TPC_main_vtx_track_nClustersPotential_TPCV1[nMaxTracks];
        int TPC_main_vtx_track_nClustersPotential_TPCV2[nMaxTracks];
        int TPC_main_vtx_track_nClustersPotential_TPCVmain[nMaxTracks];
        int TPC_main_vtx_track_nClustersPotential_TPCVgap[nMaxTracks];
        int TPC_main_vtx_track_nClustersFit_Total[nMaxTracks];
        int TPC_main_vtx_track_nClustersFit_TPCV1[nMaxTracks];
        int TPC_main_vtx_track_nClustersFit_TPCV2[nMaxTracks];
        int TPC_main_vtx_track_nClustersFit_TPCVmain[nMaxTracks];
        int TPC_main_vtx_track_nClustersFit_TPCVgap[nMaxTracks];    
        int TPC_main_vtx_track_nClustersdEdX_Total[nMaxTracks];
        int TPC_main_vtx_track_nClustersdEdX_TPCV1[nMaxTracks];
        int TPC_main_vtx_track_nClustersdEdX_TPCV2[nMaxTracks];
        int TPC_main_vtx_track_nClustersdEdX_TPCVmain[nMaxTracks];
        int TPC_main_vtx_track_nClustersdEdX_TPCVgap[nMaxTracks];
        double TPC_main_vtx_track_EnergyClusters_Total[nMaxTracks];
        double TPC_main_vtx_track_EnergyClusters_TPCV1[nMaxTracks];
        double TPC_main_vtx_track_EnergyClusters_TPCV2[nMaxTracks];
        double TPC_main_vtx_track_EnergyClusters_TPCVmain[nMaxTracks];
        double TPC_main_vtx_track_EnergyClusters_TPCVgap[nMaxTracks];
        float TPC_main_vtx_track_DCAtoVertex_X[nMaxTracks];
        float TPC_main_vtx_track_DCAtoVertex_Y[nMaxTracks];
        float TPC_main_vtx_track_DCAtoVertex_Z[nMaxTracks];
        float TPC_main_vtx_track_chi2[nMaxTracks];    
        int TPC_main_vtx_track_ndf[nMaxTracks];    
	
        int TPC_primary_vtx_track_charge[nMaxTracks];
        double TPC_primary_vtx_track_dEdx[nMaxTracks];
        int TPC_primary_vtx_track_nClusters_Total[nMaxTracks];
        int TPC_primary_vtx_track_nClusters_TPCV1[nMaxTracks];
        int TPC_primary_vtx_track_nClusters_TPCV2[nMaxTracks];
        int TPC_primary_vtx_track_nClusters_TPCVmain[nMaxTracks];
        int TPC_primary_vtx_track_nClusters_TPCVgap[nMaxTracks];
        int TPC_primary_vtx_track_nClustersPotential_Total[nMaxTracks];
        int TPC_primary_vtx_track_nClustersPotential_TPCV1[nMaxTracks];
        int TPC_primary_vtx_track_nClustersPotential_TPCV2[nMaxTracks];
        int TPC_primary_vtx_track_nClustersPotential_TPCVmain[nMaxTracks];
        int TPC_primary_vtx_track_nClustersPotential_TPCVgap[nMaxTracks];
        int TPC_primary_vtx_track_nClustersFit_Total[nMaxTracks];
        int TPC_primary_vtx_track_nClustersFit_TPCV1[nMaxTracks];
        int TPC_primary_vtx_track_nClustersFit_TPCV2[nMaxTracks];
        int TPC_primary_vtx_track_nClustersFit_TPCVmain[nMaxTracks];
        int TPC_primary_vtx_track_nClustersFit_TPCVgap[nMaxTracks];    
        int TPC_primary_vtx_track_nClustersdEdX_Total[nMaxTracks];
        int TPC_primary_vtx_track_nClustersdEdX_TPCV1[nMaxTracks];
        int TPC_primary_vtx_track_nClustersdEdX_TPCV2[nMaxTracks];
        int TPC_primary_vtx_track_nClustersdEdX_TPCVmain[nMaxTracks];
        int TPC_primary_vtx_track_nClustersdEdX_TPCVgap[nMaxTracks];
        double TPC_primary_vtx_track_EnergyClusters_Total[nMaxTracks];
        double TPC_primary_vtx_track_EnergyClusters_TPCV1[nMaxTracks];
        double TPC_primary_vtx_track_EnergyClusters_TPCV2[nMaxTracks];
        double TPC_primary_vtx_track_EnergyClusters_TPCVmain[nMaxTracks];
        double TPC_primary_vtx_track_EnergyClusters_TPCVgap[nMaxTracks];
        float TPC_primary_vtx_track_DCAtoVertex_X[nMaxTracks];
        float TPC_primary_vtx_track_DCAtoVertex_Y[nMaxTracks];
        float TPC_primary_vtx_track_DCAtoVertex_Z[nMaxTracks];
        float TPC_primary_vtx_track_chi2[nMaxTracks];    
        int TPC_primary_vtx_track_ndf[nMaxTracks];    
	
        int TPC_Multiplicity;
        int TPC_Multiplicity_all;
        int Primary_Multiplicity_all;
        int TPC_Multiplicity_Clusters_VTPC1_VTPC2;
        int TPC_Multiplicity_Clusters_All;
        float TPC_cos1;
        float TPC_sin1;
        float TPC_cos2;
        float TPC_sin2;
        
        //PSD
        int PSD_module_Number; //modules with energy
        int PSD_section_Number; //sections with energy
        int PSD_module_Number_array[nPSD_Modules]; //modules with energy
        int PSD_section_Number_array[nPSD_Modules][nPSD_Sections]; //sections with energy
        float PSD_section_slice_Energy[nPSD_Sections]; //sum over modules for a given section id
        float PSD_module_X[nPSD_Modules];
        float PSD_module_Y[nPSD_Modules];
        float PSD_module_Z[nPSD_Modules];
        float PSD_module_Energy[nPSD_Modules];
        float PSD_module_Energy_default[nPSD_Modules];
        int PSD_module_number_of_sections[nPSD_Modules];
        float PSD_section_Energy[nPSD_Modules][nPSD_Sections];
        float PSD_Energy;
        float PSD_1_Energy; //16+1
        float PSD_2_Energy; //12
        float PSD_3_Energy; //16
        
        //Main Vertex
        float Main_Vertex_X;
        float Main_Vertex_Y;
        float Main_Vertex_Z;
    
        //Primary Vertex
        float Primary_Vertex_X;
        float Primary_Vertex_Y;
        float Primary_Vertex_Z;
        bool Has_Primary_Vertex;

        //Additionary Vertex

        float pVtx_Z;
    
    //Triggers    
    unsigned int BPD_Status[nRawReco][nBPD][nBPDcomponents];
    float BPD_Position[nRawReco][nBPD][nBPDcomponents];
    float BPD_PositionError[nRawReco][nBPD][nBPDcomponents];
    float BPD_Z[nRawReco][nBPD][nBPDcomponents];
    float BPD_RMS[nRawReco][nBPD][nBPDcomponents];
    float BPD_Maximum[nRawReco][nBPD][nBPDcomponents];
    float BPD_Charge[nRawReco][nBPD][nBPDcomponents];
    float BPD_SumOfAll[nRawReco][nBPD][nBPDcomponents];
    float triggersADC[nRawReco][nTriggers_Simple];
    bool isTriggers_Simple[nRawReco][nTriggers_Simple];
    bool isTriggers_Combined[nRawReco][nTriggers_Combined];
    
    //Beam
    float Beam_Momentum[nRawReco][nBeamComponents];
    float Beam_Fitted2DLineXZ[nRawReco][nBeamComponents];
    float Beam_Fitted2DLineYZ[nRawReco][nBeamComponents];
    int Beam_Status[nRawReco];
//     float beamPosition[nRawReco][nBeamComponents];

    //WFA
    float WFA_TimeStructure[nTriggers_Simple][nMaxWFAsignals];
    unsigned int WFA_NumberOfSignalHits[nTriggers_Simple];
    
    //FitVertex (primary vertex)
    float FitVertexX;
    float FitVertexY;
    float FitVertexZ;
//     EFitQuality FitVertexQ;
    int FitVertexQ;
    int Primary_Vertex_Q;
//     		Vertex* pPrimaryVertex = &event.GetRecEvent().GetPrimaryVertex(rec::VertexConst::ePrimaryFitZ);
// 		const utl::Point& Vertex = pPrimaryVertex->GetPosition();
// 		myFitVtxZ = Vertex.GetZ();
//         myFitVtxX = Vertex.GetX();
//         myFitVtxY = Vertex.GetY();
// 	myFitVtxQ = Vertex.GetFitQuality

    
    
//	TODO ask:
//     	??
//     	WFA, 4μs, WFA interaction, 25μs
//     	ntf/nto>0.25?
//     	Fit quality?
// 	//TODO analysis:
//  	Vertex track status+has track
// 	I Min number of clusters in VTPCs - 15
// 	I Min number of clusters in all TPCs - 30.
// 	I Impact parameter, X < 4cm, Y < 2cm
    //Old style variables for triggers
    bool T1;
    bool T2;
    bool T3;
    bool T4;
    float S1;
    float S2;
    float S3;
    float S;
    float V1;
    
    TFile *outputFile;

    REGISTER_MODULE("FlowQA", FlowQA,"$Id: FlowQA.h $");
  };

}


#endif
