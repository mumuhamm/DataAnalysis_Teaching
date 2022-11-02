void example(){

	//Simple functions to test 3d plotting capability
  TF3 *func1 = new TF3("func1","[0]*exp(-0.5*((x-[1])/[2])**2)*[0]*exp(-0.5*((y-[1])/[2])**2)*[0]*exp(-0.5*((z-[1])/[2])**2)");
  TF3 *func2 = new TF3("func2","1/(x^2+y^2+z^2)^(1/2)");
 

	//Spherical Harmonics:  Ylm
  TF3 *Y00 = new TF3("Y00","1/(x^2+y^2+z^2)");
  TF3 *Y10 = new TF3("Y10","1/(x^2+y^2+z^2)*((3/(4*3.14159))*(cos(atan((x^2+y^2)^(1/2)/z)))^2)");
  TF3 *Y11 = new TF3("Y11","1/(x^2+y^2+z^2)*((3/(8*3.14159))^(1/2)*(sin(atan((x^2+y^2)^(1/2)/z)))^2*(cos(atan(x/y)))^2)");
  TF3 *Y20 = new TF3("Y20","1/(x^2+y^2+z^2)*((5/(16*3.14159))^(1/2)*(3*(cos(atan((x^2+y^2)^(1/2)/z)))^2-1))^2");
  TF3 *Y21 = new TF3("Y21","1/(x^2+y^2+z^2)*((15/(8*3.14159))^(1/2)*(sin(atan((x^2+y^2)^(1/2)/z))*cos(atan((x^2+y^2)^(1/2)/z))*cos(acos(x/(x^2+y^2)^(1/2)))))^2");
  TF3 *Y22 = new TF3("Y22","1/(x^2+y^2+z^2)*((15/(32*3.14159))^(1/2)*((sin(atan((x^2+y^2)^(1/2)/z)))^2*(cos(acos(x/(x^2+y^2)^(1/2))))^2))^2");

	//Wave functions for 3D Simple Harmonic Oscillator:  SHO_nlm
 TF3 *SHO_000 = new TF3("SHO_000","((2*[0]^(3/2)/3.14159^(1/4))*exp(-[0]^2*(x^2+y^2+z^2)/2))^2");
 TF3 *SHO_110 = new TF3("SHO_110","((2*[0]^(3/2)*(2/3)^(1/2)/3.14159^(1/4))*([0]*(x^2+y^2+z^2)^(1/2))*exp(-[0]^2*(x^2+y^2+z^2)/2))^2*(((3/(4*3.14159))*(cos(atan((x^2+y^2)^(1/2)/z)))^2))");
 TF3 *SHO_111 = new TF3("SHO_111","((2*[0]^(3/2)*(2/3)^(1/2)/3.14159^(1/4))*([0]*(x^2+y^2+z^2)^(1/2))*exp(-[0]^2*(x^2+y^2+z^2)/2))^2*(((3/(8*3.14159))^(1/2)*(sin(atan((x^2+y^2)^(1/2)/z)))^2*(cos(atan(x/y)))^2))");
 TF3 *SHO_200 = new TF3("SHO_200","((2*[0]^(3/2)*(2/3)^(1/2)/3.14159^(1/4))*(3/2 -[0]^2*(x^2+y^2+z^2))*exp(-[0]^2*(x^2+y^2+z^2)/2))^2");
 TF3 *SHO_220 = new TF3("SHO_220","((4*[0]^(3/2)*(1/15)^(1/2)/3.14159^(1/4))*([0]^2*(x^2+y^2+z^2))*exp(-[0]^2*(x^2+y^2+z^2)/2))^2*((5/(16*3.14159))^(1/2)*(3*(cos(atan((x^2+y^2)^(1/2)/z)))^2-1))^2");
 TF3 *SHO_221 = new TF3("SHO_221","((4*[0]^(3/2)*(1/15)^(1/2)/3.14159^(1/4))*([0]^2*(x^2+y^2+z^2))*exp(-[0]^2*(x^2+y^2+z^2)/2))^2*((15/(8*3.14159))^(1/2)*(sin(atan((x^2+y^2)^(1/2)/z))*cos(atan((x^2+y^2)^(1/2)/z))*cos(acos(x/(x^2+y^2)^(1/2)))))^2");
 TF3 *SHO_222 = new TF3("SHO_222","((4*[0]^(3/2)*(1/15)^(1/2)/3.14159^(1/4))*([0]^2*(x^2+y^2+z^2))*exp(-[0]^2*(x^2+y^2+z^2)/2))^2*((15/(32*3.14159))^(1/2)*((sin(atan((x^2+y^2)^(1/2)/z)))^2*(cos(2*acos(x/(x^2+y^2)^(1/2))))))^2");

	//Set Parameters for functions to be plotted
 SHO_000->SetParameter(0,3);
 SHO_110->SetParameter(0,3);
 SHO_111->SetParameter(0,3);
 SHO_200->SetParameter(0,3);
 SHO_220->SetParameter(0,3);
 SHO_221->SetParameter(0,3);
 SHO_222->SetParameter(0,3);
 
	//3D Histogram
 //  TH3F *histo = new TH3F("histo","3d hist",100,-1,1,100,-1,1,100,-1,1);
  TH3F *histo = new TH3F("histo","3d hist",50,-1,1,50,-1,1,50,-1,1);

	//Fill the histogram randomly according to the distribution
	//      defined by the indicated function.
	
    //Change the expression in quotes to switch the function plotted
	//Change the number to alter the number of points plotted

	histo->FillRandom("SHO_222",100000);
	// histo->FillRandom("Y10",100000);

 //Ttree method


 struct mystruct{
   Double_t xval;
   Double_t yval;
   Double_t zval;
   Int_t weight;
 };
 mystruct stuff;


 TTree *tree = new TTree("tree","my tree");
 tree->Branch("mybranch",&stuff.xval,"xval/D:yval:zval:weight/I");



 for(Int_t i=0;i<histo->GetXaxis()->GetNbins();i++){

for(Int_t j=0;j<histo->GetYaxis()->GetNbins();j++){
for(Int_t k=0;k<histo->GetZaxis()->GetNbins();k++){
    if(histo->GetBinContent(i,j,k)>0){

      //  stuff.xval=histo->GetXaxis()->GetBinCenter(i);
      // stuff.yval=histo->GetYaxis()->GetBinCenter(j);
      // stuff.zval=histo->GetZaxis()->GetBinCenter(k);
   Double_t centerx=histo->GetXaxis()->GetBinCenter(i);
   Double_t centery=histo->GetYaxis()->GetBinCenter(j);
   Double_t centerz=histo->GetZaxis()->GetBinCenter(k);
    Double_t widthx=histo->GetXaxis()->GetBinWidth(i);
   Double_t widthy=histo->GetYaxis()->GetBinWidth(j);
   Double_t widthz=histo->GetZaxis()->GetBinWidth(k);
 
   // if(histo->GetBinContent(i,j,k)>0){stuff.weight = histo->GetBinContent(i,j,k);}else{stuff.weight =0;}
    for(Int_t v =0;v<10;v++){

      stuff.xval = gRandom->Gaus(centerx,widthx);
      stuff.yval = gRandom->Gaus(centery,widthy);
      stuff.zval = gRandom->Gaus(centerz,widthz);

    stuff.weight = histo->GetBinContent(i,j,k);
       tree->Fill();
   }
   }
  //stuff.weight = histo->GetBinContent(i,j,k);
 
   // tree->Fill();
  // }
 }
 }


}


 TCanvas *c1 = new TCanvas("c1","c1");
 c1->cd();
 tree->Print();
 tree->Draw("xval:yval:zval:weight","");

}
