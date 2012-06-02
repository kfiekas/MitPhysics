#include "graph_2catVBF.C"
#include "graph_3catVBF.C"
#include "graph_4catVBF.C" 
#include "graph_5catVBF.C"


void makePlot() {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.035, "XYZ");

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleOffset(1., "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");

  gStyle->SetPalette(1);

  //minimum limit line for each case of number of categories
  TH1D* _1Cat = new TH1D("_1Cat","",1,-0.5,1.);
  _1Cat->SetBinContent(1,1.62);
  _1Cat->SetLineStyle(kDashed);
  _1Cat->SetLineColor(kBlack);
  _1Cat->SetMinimum(0.95);
  _1Cat->SetMaximum(1.7);
  _1Cat->GetXaxis()->SetTitle("Category Split (BDT)");
  _1Cat->GetYaxis()->SetTitle("#sigma_{95%}/#sigma_{SM}");

  TH1D* _2Cat = new TH1D("_2Cat","",1,-0.5,1.);
  _2Cat->SetLineStyle(kDashed);
  _2Cat->SetLineColor(kBlue);
  _2Cat->SetBinContent(1,1.1554);

  TH1D* _3Cat = new TH1D("_3Cat","",1,-0.5,1.);
  _3Cat->SetLineStyle(kDashed);
  _3Cat->SetLineColor(kGreen+2);
  _3Cat->SetBinContent(1,1.0675);

  TH1D* _4Cat = new TH1D("_4Cat","",1,-0.5,1.);
  _4Cat->SetLineStyle(kDashed);
  _4Cat->SetLineColor(kMagenta);
  _4Cat->SetBinContent(1,1.0358);

  TH1D* _5Cat = new TH1D("_5Cat","",1,-0.5,1.);
  _5Cat->SetLineStyle(kDashed);
  _5Cat->SetLineColor(kRed);
  _5Cat->SetBinContent(1,1.0187);

  /* TH1D* _6Cat = new TH1D("_6Cat","",1,-0.5,1.);
  _6Cat->SetLineStyle(kDashed);
  _6Cat->SetLineColor(kRed);
  _6Cat->SetBinContent(1,);*/ 
  
  /*TH1D* lp4Cat = new TH1D("lp4Cat","",1,.-0.7,1.);
    lp4Cat->SetBinContent(1,2.1741);
    lp4Cat->SetLineColor(kGreen+2);
    lp4Cat->SetLineStyle(kSolid);
    lp4Cat->SetLineWidth(2.);*/
  
  
  // IMMPORT RESULT HERE
  // limit VS newly added cut for each case of number of categories
  TGraph* obsCat2      = NULL;
  graph_2catVBF(obsCat2);
  obsCat2->SetLineColor(kBlue);
  obsCat2->SetLineStyle(kSolid);

  TGraph* obsCat3      = NULL;
  graph_3catVBF(obsCat3);
  obsCat3->SetLineColor(kGreen+2);
  obsCat3->SetLineStyle(kSolid);

  TGraph* obsCat4      = NULL;
  graph_4catVBF(obsCat4);
  obsCat4->SetLineColor(kMagenta);
  obsCat4->SetLineStyle(kSolid);

  TGraph* obsCat5     = NULL;
  graph_5catVBF(obsCat5);
  obsCat5->SetLineColor(kRed);
  obsCat5->SetLineStyle(kSolid);
  obsCat5->SetLineWidth(1.5);

  /*TGraph* obsCat6      = NULL;
  graph_6catNew(obsCat6);
  obsCat6->SetLineColor(kBlue);
  obsCat6->SetLineStyle(kSolid);
  obsCat6->SetLineWidth(1.5);*/
  
  // Catrgpry lines: pependicular line for the minimum limit for each case of number of categories
  double x2[2];
  double y2[2];
  x2[0]=0.545;
  x2[1]=0.545;
  y2[0]=-10;
  y2[1]=1.1554;
  TGraph* cat2L = new TGraph(2,x2,y2);

  double x3[2];
  double y3[2];
  x3[0]=0.89;
  x3[1]=0.89;
  y3[0]=-10;
  y3[1]=1.0675;
  TGraph* cat3L = new TGraph(2,x3,y3);

  double x4[2];
  double y4[2];
  x4[0]=0.05;
  x4[1]=0.05;
  y4[0]=-10.;
  y4[1]=1.0358;
  TGraph* cat4L = new TGraph(2,x4,y4);

  double x5[2];
  double y5[2];
  x5[0]=0.74;
  x5[1]=0.74;
  y5[0]=-10.;
  y5[1]=1.0187;
  TGraph* cat5L = new TGraph(2,x5,y5);

  //double x6[2];
  //double y6[2];
  //x6[0]=0.75625;
  //x6[1]=0.75625;
  //y6[0]=-10.;
  //y6[1]=1.0199;
  //TGraph* cat5L = new TGraph(2,x6,y6);
 
  cat2L->SetLineColor(kBlue);
  cat3L->SetLineColor(kGreen+2);
  cat4L->SetLineColor(kMagenta);
  cat5L->SetLineColor(kRed);
  //cat6L->SetLineColor(kRed);
 

  cat2L->SetLineWidth(2);
  cat3L->SetLineWidth(2);
  cat4L->SetLineWidth(2);
  cat5L->SetLineWidth(2);
  //cat6L->SetLineWidth(2);

  _1Cat->Draw();
  _2Cat->Draw("SAME");
  _3Cat->Draw("SAME");
  _4Cat->Draw("SAME");
  _5Cat->Draw("SAME");
  //_6Cat->Draw("SAME");


  cat2L->Draw("LP");
  cat3L->Draw("LP");
  cat4L->Draw("LP");
  cat5L->Draw("LP");
  //cat6L->Draw("LP");

  //lp4Cat->Draw("SAME");

  obsCat2->Draw("LP");
  obsCat3->Draw("LP");
  obsCat4->Draw("LP");
  obsCat5->Draw("LP");
  //obsCat6->Draw("LP");

  /*TLatex* text_LP=new TLatex(3.5,23.5,"2#eta x 2R9 Cats");
  text_LP->SetNDC();
  text_LP->SetTextAlign(13);
  text_LP->SetX(0.134);//(0.940);
  text_LP->SetY(0.564);
  text_LP->SetTextFont(42);
  text_LP->SetTextColor(kGreen+2);
  text_LP->SetTextSize(0.03);// dflt=28

  text_LP->Draw();*/

  TLatex* text_1Cat=new TLatex(3.5,23.5,"1 MVA Cat");
  text_1Cat->SetNDC();
  text_1Cat->SetTextAlign(13);
  text_1Cat->SetX(0.134);//(0.940);
  text_1Cat->SetY(0.85);
  text_1Cat->SetTextFont(40);
  text_1Cat->SetTextSize(0.028);// dflt=28
  text_1Cat->SetTextColor(kBlack);
  text_1Cat->Draw();

  TLatex* text_2Cat=new TLatex(3.5,23.5,"2 MVA Cats");
  text_2Cat->SetNDC();
  text_2Cat->SetTextAlign(13);
  text_2Cat->SetX(0.134);//(0.940);
  text_2Cat->SetY(0.35);
  text_2Cat->SetTextFont(40);
  text_2Cat->SetTextSize(0.028);// dflt=28
  text_2Cat->SetTextColor(kBlue);
  text_2Cat->Draw();

  TLatex* text_3Cat=new TLatex(3.5,23.5,"3 MVA Cats");
  text_3Cat->SetNDC();
  text_3Cat->SetTextAlign(13);
  text_3Cat->SetX(0.134);//(0.940);
  text_3Cat->SetY(0.26);
  text_3Cat->SetTextFont(40);
  text_3Cat->SetTextSize(0.028);// dflt=28
  text_3Cat->SetTextColor(kGreen+2);
  text_3Cat->Draw();

  TLatex* text_4Cat=new TLatex(3.5,23.5,"4 MVA Cats");
  text_4Cat->SetNDC();
  text_4Cat->SetTextAlign(13);
  text_4Cat->SetX(0.134);//(0.940);
  text_4Cat->SetY(0.22);
  text_4Cat->SetTextFont(40);
  text_4Cat->SetTextSize(0.028);// dflt=28
  text_4Cat->SetTextColor(kMagenta);
  text_4Cat->Draw();

  TLatex* text_5Cat=new TLatex(3.5,23.5,"5 MVA Cats");
  text_5Cat->SetNDC();
  text_5Cat->SetTextAlign(13);
  text_5Cat->SetX(0.134);//(0.940);
  text_5Cat->SetY(0.16);
  text_5Cat->SetTextFont(40);
  text_5Cat->SetTextSize(0.028);//dflt=28
  text_5Cat->SetTextColor(kRed);
  text_5Cat->Draw();

  TLatex* text_6Cat=new TLatex(3.5,23.5,"5 MVA Cats");
  text_6Cat->SetNDC();
  text_6Cat->SetTextAlign(13);
  text_6Cat->SetX(0.134);//(0.940);
  text_6Cat->SetY(0.16);
  text_6Cat->SetTextFont(42);
  text_6Cat->SetTextSize(0.028);//dflt=28
  text_6Cat->SetTextColor(kRed);
  text_6Cat->Draw();



  /*TLatex* text_Scan2=new TLatex(3.5,23.5,"Scan 1#rightarrow2 Cats");
  text_Scan2->SetNDC();
  text_Scan2->SetTextAlign(13);
  text_Scan2->SetX(0.414);//(0.940);
  text_Scan2->SetY(0.648);
  text_Scan2->SetTextFont(42);
  text_Scan2->SetTextColor(kBlue);
  text_Scan2->SetTextSize(0.03);// dflt=28

  //text_Scan2->Draw();

  TLatex* text_Scan3=new TLatex(3.5,23.5,"Scan 2#rightarrow3 Cats");
  text_Scan3->SetNDC();
  text_Scan3->SetTextAlign(13);
  text_Scan3->SetX(0.517);//(0.940);
  text_Scan3->SetY(0.383);
  text_Scan3->SetTextFont(42);
  text_Scan3->SetTextColor(kBlue);
  text_Scan3->SetTextSize(0.03);// dflt=28

  //text_Scan3->Draw();

  TLatex* text_Scan4=new TLatex(3.5,23.5,"Scan 3#rightarrow4 Cats");
  text_Scan4->SetNDC();
  text_Scan4->SetTextAlign(13);
  text_Scan4->SetX(0.549);//(0.940);
  text_Scan4->SetY(0.286);
  text_Scan4->SetTextFont(42);
  text_Scan4->SetTextColor(kBlue);
  text_Scan4->SetTextSize(0.03);// dflt=28

  //text_Scan4->Draw();

  TLatex* text_Scan5=new TLatex(3.5,23.5,"Scan 3#rightarrow4 Cats");
  text_Scan5->SetNDC();
  text_Scan5->SetTextAlign(13);
  text_Scan5->SetX(0.549);//(0.940);
  text_Scan5->SetY(0.286);
  text_Scan5->SetTextFont(42);
  text_Scan5->SetTextColor(kBlue);
  text_Scan5->SetTextSize(0.03);// dflt=28*/

  //text_Scan5->Draw();


  /*TLatex* text_Split1=new TLatex(3.5,23.5,"1#rightarrow2 (0.730)");
  text_Split1->SetNDC();
  text_Split1->SetTextAlign(13);
  text_Split1->SetX(0.747);//(0.940);
  text_Split1->SetY(0.118);
  text_Split1->SetTextFont(42);
  text_Split1->SetTextAngle(90);
  text_Split1->SetTextColor(kRed);
  text_Split1->SetTextSize(0.025);// dflt=28

  //text_Split1->Draw();

  TLatex* text_Split21=new TLatex(3.5,23.5,"2#rightarrow3 (0.875)");
  text_Split21->SetNDC();
  text_Split21->SetTextAlign(13);
  text_Split21->SetX(0.819);//(0.940);
  text_Split21->SetY(0.118);
  text_Split21->SetTextFont(42);
  text_Split21->SetTextAngle(90);
  text_Split21->SetTextColor(kRed);
  text_Split21->SetTextSize(0.025);// dflt=28

  //text_Split21->Draw();

  TLatex* text_Split22=new TLatex(3.5,23.5,"3#rightarrow4 (0.450)");
  text_Split22->SetNDC();
  text_Split22->SetTextAlign(13);
  text_Split22->SetX(0.619);//(0.940);
  text_Split22->SetY(0.118);
  text_Split22->SetTextFont(42);
  text_Split22->SetTextAngle(90);
  text_Split22->SetTextColor(kRed);
  text_Split22->SetTextSize(0.025);// dflt=28

  //text_Split22->Draw();*/

  TLatex* text=new TLatex(3.5,23.5,"CMS Preliminary");
  text->SetNDC();
  //text->SetTextAlign(13);
  text->SetX(0.1);//(0.940);
  text->SetY(0.92);
  text->SetTextFont(42);
  text->SetTextSize(0.035);// dflt=28

  TLatex*  text2 = new TLatex(3.5,23.5,"2012Jan16 MC (L = 5.089 fb^{-1})");
  text2->SetNDC();
  //text2->SetTextAlign(33);
  text2->SetX(0.4);//(0.940);
  text2->SetY(0.92);
  text2->SetTextFont(42);
  text2->SetTextSize(0.035);// dflt=28

  TLatex* text3=new TLatex(3.5,23.5,"M_{H}=125 GeV");
  text3->SetNDC();
  //text3->SetTextAlign(33);
  text3->SetX(0.76);//(0.940);
  text3->SetY(0.92);
  text3->SetTextFont(42);
  text3->SetTextSize(0.035);// dflt=28

  text->Draw();
  text2->Draw();
  text3->Draw();

  TLegend* leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);  

  leg->AddEntry(obsCat2," 2 MVA Cats ","L");
  leg->AddEntry(obsCat3," 3 MVA Cats ","L");
  leg->AddEntry(obsCat4," 4 MVA Cats ","L");
  leg->AddEntry(obsCat5," 5 MVA Cats ","L");
  //leg->AddEntry(obsCat6," 6 MVA Cats ","L");

  leg->Draw();

  return;
}
