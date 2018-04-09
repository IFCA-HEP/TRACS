/*
 * @ Copyright 2014-2017 CERN and Instituto de Fisica de Cantabria - Universidad de Cantabria. All rigths not expressly granted are reserved [tracs.ssd@cern.ch]
 * This file is part of TRACS.
 *
 * TRACS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the Licence.
 *
 * TRACS is distributed in the hope that it will be useful , but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with TRACS. If not, see <http://www.gnu.org/licenses/>
 */

// In general this is a macro to convolute 2 histograms, quite general
// The best is to use it as an executable, the histos are saved into a root file
// 
//   H1DConvolution [Cend in pF] [Pulse duration in ns]
//
// Using it as a script:
//
// Example: Convolution of a TCT-like signal with the transfer function of the amplifier
// root -l $SRC/Centered_100ps_TransferFunction_Cividec_06052014.root
// TH1D *tf=_file0->Get("shtf")
// TH1D *tct=_file0->Get("TCT signal")
// .L $SRC/H1DConvolution.C++ 
// TH1D *hconv=H1DConvolution(tf,tct)

#include <cstdarg>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <mutex>          // std::mutex


#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TAttLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "THStack.h"
#include "TF1.h"


#define EXE 1      //1 for shared lib inside root. Note that still a ".o" can be produced in this mode
                   //by compiling the source, as needed for plotvd.C
                   
                   //2 for executable outside root

#define CPF 2.0		   

#define CTRLPLOT 0      //1 if intermediate control plots wanted

//See "Visual explanations of convolution" in http://en.wikipedia.org/wiki/Convolution
TH1D *H1DConvolution( TH1D *htf , TH1D *htct , Double_t Cend=0. , int tid=0) ; 
TH1D *H1DConvolution( TH1D *htct  , Double_t Cend=0. , int tid=0) ; 
TH1D *LPFilter( TH1D *htf , Double_t Cend ) ; 

std::mutex mtx_conv;           // mutex for critical section
int count;

TH1D *LPFilter( TH1D *hin , Double_t Cend  ) {

    if (Cend==-1.0) Cend = CPF ; 
    Double_t RCns = 50. * Cend * 1.e-3 ;
    Double_t At = hin->GetBinCenter(2) - hin->GetBinCenter(1) ;
    Double_t alfa = At/(RCns + At) ;
    TH1D *hout = (TH1D *) hin->Clone(); hout->Reset();
    hout->SetTitle( Form("C=%d pF",TMath::Nint(Cend)) ) ;
    hout->SetName( Form("C=%d pF",TMath::Nint(Cend)) ) ;

    hout->SetBinContent(1, hin->GetBinContent(1) );
    for ( Int_t i=2 ; i<=hin->GetNbinsX() ; i++ ) {
      Double_t val = (1.0-alfa)*hout->GetBinContent(i-1) + alfa*hin->GetBinContent(i) ; 
      hout->SetBinContent(i,val);
    }
    
    return hout ;
    
}

TH1D *H1DConvolution( TH1D *htf , TH1D *htct , Double_t Cend , int tid) { 
      
   //------------>Both input histograms should have the same bin width<-----------------
   
   //Here you can apply an extra LPFiltering
   //if (Cend!=0) htct = LPFilter( htct , Cend ); 

   //Convolute (commutative)
   //C(t) = Int[ tct(x) transferfunction(t-x) dx ]
   Double_t bw = htct->GetBinCenter(2) - htct->GetBinCenter(1);
   //Double_t bw = htct->GetXaxis()->GetBinCenter(2) - htct->GetXaxis()->GetBinCenter(1);
   //TF1 *f1 = new TF1("f1","abs(sin(x)/x)*sqrt(x)",0,10);
   //   float r = f1->GetRandom();
    TString tftit, tfname;
    tftit.Form("hConv_%d_%d", tid, count);
    tfname.Form("conv_%d_%d", tid, count);
   TH1D *hConv = new TH1D(tftit,tfname,2*htct->GetNbinsX(),-htct->GetNbinsX()*bw,htct->GetNbinsX()*bw);
   
   
   //The convoluted response to the TCT signal is going to be another histogram sized similar to htct
   Int_t Ntf = htf->GetNbinsX() , Ntct = htct->GetNbinsX();
   
   //Create the reverse histogram of the transfer function
   TH1D *hinv = (TH1D *) htf->Clone(); hinv->Reset();
   for (Int_t j=1; j<= Ntf ; j++) hinv->SetBinContent( j , htf->GetBinContent(Ntf-j+1) );  
  // std::basic_string<char> t= std::to_string(tid);
   TH1D *hg = (TH1D *) htct->Clone(); 
   
   hinv->Draw();
   #if CTRLPLOT == 1
   gPad->Print( "hgcontrol.pdf[" ) ; gPad->Print( "hgcontrol.pdf" );
   #endif
   
   for ( Int_t i=1 ; i<=2*Ntct ; i++ ) { 
          
     //Create a shifted histogram version of the inverse
     hg->Reset();
     //for (Int_t j=TMath::Nint(-0.5*Ntf); j<= TMath::Nint(0.5*Ntf) ; j++) hg->SetBinContent( i-j , hinv->GetBinContent( j+TMath::Nint(0.5*Ntf) ) );  
     if ( i<=Ntf ) {
       //Histogram is shifting in from the left
       for (Int_t j=1; j<=i ; j++) hg->SetBinContent( j , hinv->GetBinContent( Ntf-i+j ) );
     } else {
       //Histogram is shifting out. Leaving from the right
       Int_t cont=1 ;
       for (Int_t j=i-Ntf+1; j<=2*Ntf ; j++) {
         hg->SetBinContent( j , hinv->GetBinContent( cont ) );
	 cont++;
       }
  
     }
     
     
     //Multiply f(tau)*g(t-tau)
     hg->Multiply( htct );
     
     //Double_t fxg = hg->Integral("width");
     Double_t fxg = hg->Integral();
     
     hConv->SetBinContent(i,fxg);
     
     #if CTRLPLOT==1
       THStack *hst=new THStack("hst","conv");
       hConv->SetLineColor(2);hConv->SetLineWidth(2);
       hst->Add(hg) ; hst->Add(hConv);
       hst->Draw("nostack");
       if   (i==2*Ntct) { 
         gPad->Print( "hgcontrol.pdf" )  ; gPad->Print( "hgcontrol.pdf]" ) ; 
       } else if (i%10==0)  gPad->Print( "hgcontrol.pdf" )  ;
     #endif

   }

   gStyle->SetOptStat(0);
   gStyle->SetHistLineWidth(2);
   //hConv->SetLineColor(kRed) ;htct->SetLineColor(kBlack) ;htf->SetLineColor(kBlue) ;
   
   THStack *hs = new THStack();
   //hs->Add(htf);
   hs->Add(htct);
   hs->Add(hConv);
   TCanvas *c1=new TCanvas(); c1->cd();
   hs->Draw("nostack");
   //hs->GetXaxis()->SetRangeUser(-2.,10.);
   hs->GetXaxis()->SetTitle("Time [ns]") ;

   c1->SetGrid(1);

   TLegend* legend = c1->BuildLegend();
   legend->Draw();
   c1->Update();
   //   NOPDF for the moment
//   c1->Print( "convolution.pdf" );
 //too many files...
/*
    tftit.Form("conv_%d_%d.root", tid, count);
   TFile *f=new TFile(tftit,"UPDATE");
   hConv->Write();
   //f->Close();
   delete f;

   #if EXE==1
      tftit.Form("convolution_%d_%d.root", tid, count);
      TFile *fout=new TFile(tftit,"UPDATE");
      htct->Write();
      hConv->Write();
      //fout->Close();
      delete fout;
   #endif
  count++; //for naming purposes
*/

   return hConv;  

}

TH1D *H1DConvolution( TH1D *htct , Double_t Cend, int tid) { 
   
   //mtx_conv.lock();
   TFile  *ftf = new TFile( "Centered_100ps_TransferFunction_Cividec_06052014.root");
   TH1D   *htf = (TH1D *) ftf->Get("shtf");
   TH1D *hConv = H1DConvolution(htf,htct,Cend,tid);
   //mtx_conv.unlock();
   return hConv;
   
}


#if EXE==2
int main (int argc , char *argv[] ) {

   Double_t Cend = 0. , PulseDuration = 0.5 ;
   if (argc==2) Cend = atof( argv[1] );
   if (argc==3) {
     Cend          = atof( argv[1] );
     PulseDuration = atof( argv[2] );
   }
   
   //Open histograms: Transfer Function
   //TFile *ftf=new TFile("/home/mfg/hpk/data/PulseGenerator/100ps_TransferFunction_Cividec_06052014.root");
   //TFile *ftf=new TFile("/home/mfg/hpk/ETCT-Analyse/cpp/Miteq_NetworkAnalyzer.root");
   TFile *ftf=new TFile("~/etct/ETCT-Analyse/mfg/Centered_100ps_TransferFunction_Cividec_06052014.root");
   TKey *key ; 
   TIter nextkey( ftf->GetListOfKeys() );
   key = (TKey*)nextkey() ;

   TObject *obj = key->ReadObj();
   TH1D *htf = (TH1D*) obj;
   htf->SetTitle("Transfer function");
   
   //Center TF, by padding it with zeroes
   Double_t Xmin = htf->GetXaxis()->GetXmin() , Xmax = htf->GetXaxis()->GetXmax() , nXmin, nXmax;
   if ( TMath::Abs(Xmin) > TMath::Abs(Xmax) ) { nXmin = Xmin                 ; nXmax = TMath::Abs(Xmin) ; }
   if ( TMath::Abs(Xmin) < TMath::Abs(Xmax) ) { nXmin = -1.*TMath::Abs(Xmax) ; nXmax = Xmax             ; } 
   Int_t Nbins= TMath::Nint( (nXmax-nXmin)/(htf->GetBinCenter(2)-htf->GetBinCenter(1)) );
   
   TH1D *shtf = new TH1D( "shtf" , "Trans. Funct." , Nbins , nXmin , nXmax );
   
   for (Int_t i=1;i<=Nbins;i++) {
     if ( shtf->GetBinCenter(i)>=Xmin && shtf->GetBinCenter(i)<=Xmax  ) shtf->SetBinContent( i , htf->GetBinContent( htf->FindBin(shtf->GetBinCenter(i)) )) ;
   }   
   
      
   //Fake signal
   TH1D *htct=(TH1D *) shtf->Clone("TCT signal"); htct->Reset() ; htct->SetTitle("Signal");  
   
   for (Int_t i=1;i<=Nbins;i++) {

     //box-like pulse (signal coming from HVCMOS)
     if   ( htct->GetBinCenter(i)>=0. && htct->GetBinCenter(i)<=PulseDuration ) htct->SetBinContent( i , 0.05 ) ;
     
     //TCT-like pulse (invented)
     //if   ( htct->GetBinCenter(i)>=0. && htct->GetBinCenter(i)<=0.5 ) htct->SetBinContent( i , 0.05-5.6e-3*htct->GetBinCenter(i) ) ;
     //else htct->SetBinContent( i , 0.0 ) ;
   }
   
   TH1D *htctLP = LPFilter( htct , Cend );
   
   TFile *fout=new TFile( Form("Convolution_%dpF_%dns.root",TMath::Nint(Cend),TMath::Nint(PulseDuration)),"RECREATE" );
   //TFile *fout=new TFile("Centered_signaltctlike_100ps_TransferFunction_Cividec_06052014.root","RECREATE");
   shtf->Write();
   htct->Write();
   htctLP->Write();

   TH1D *hconv = H1DConvolution( shtf, htctLP , 0);
   
   hconv->Write();
   fout->Close();
   
   #ifdef FAKEINPUT
     //Fake signal
     TH1D *hb = new TH1D("hb","hb"   , 100 , -6 , 6);

     //Fill the box
     for (Int_t i=1;i<=100;i++) {
       hb->SetBinContent( i , 0.0 ) ;
       if ( i>42&&i<58 ) hb->SetBinContent( i , 1.0 ) ;
     }

     //Gaussian as transfer function
     TH1D *hg = new TH1D("hg","Transf. function"   , 100 , -6 , 6);   
     for (Int_t i=1;i<=100;i++) hg->SetBinContent( i , TMath::Gaus(-6+(i-1)*12./100.,0.,.5,kTRUE) );

     H1DConvolution( hg, hb ,0);
   #endif
   
   
   
   return 0 ;
}
#endif
