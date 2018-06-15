#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TH3D.h"

#include <iostream>
#include <map>

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"

using namespace std;

int main(int argc, char* argv[])
{
    const char* inputroot = argv[1];
    int cluster_check=-1; // specify cluster_id to check 
    if(argc==3) cluster_check=atoi(argv[2]);

    // check
    TH3D* hrec = new TH3D("hrec","",25,0,256,25,-115,117,100,0,1037);
    TPolyLine3D *leval1 = new TPolyLine3D(2);
    TPolyLine3D *leval2 = new TPolyLine3D(2);

    // read in the point cloud for each cluster_id
    TFile* f = new TFile(inputroot, "READ");
    TTree* clusters = (TTree*)f->Get("T_cluster");
    Double_t x=0;
    Double_t y=0;
    Double_t z=0;
    Int_t cluster_id=0;

    clusters->SetBranchAddress("x",&x);
    clusters->SetBranchAddress("y",&y);
    clusters->SetBranchAddress("z",&z);
    clusters->SetBranchAddress("cluster_id",&cluster_id);

    std::map<int, std::vector<double>> xpt;
    std::map<int, std::vector<double>> ypt;
    std::map<int, std::vector<double>> zpt;

    for(int i=0; i<clusters->GetEntries();i++)
    {
        clusters->GetEntry(i);
        auto it = xpt.find(cluster_id);
        if( it == xpt.end() ){
            std::vector<double> xvec;
            xvec.push_back(x);
            std::vector<double> yvec;
            yvec.push_back(y);
            std::vector<double> zvec;
            zvec.push_back(z);

            xpt[cluster_id] = xvec;
            ypt[cluster_id] = yvec;
            zpt[cluster_id] = zvec;
        }
        else{
            xpt[cluster_id].push_back(x);
            ypt[cluster_id].push_back(y);
            zpt[cluster_id].push_back(z);
        }
    }
    //cout<<clusters->GetEntries()<<endl;

    // PCA (Priciple Component Analysis): Matrix Decomposition --> Coordinate Rotation
    // n-dimentional points in nxn position covariance matrix;
    // diagonalization (rotation);
    // the biggest eigen value corresponds to the main axis --> primary direction of the point cloud
    
    // Step 1: calculate average position
    // Step 2: construct covariance matrix
    // merge the two steps into a single loop of the point cloud

    //cout<<xpt.size()<<endl;
    for(int i=0; i<xpt.size(); i++)
    { // each cluster
        double npts = xpt[i].size();
        double xx=0;
        double yy=0;
        double zz=0;
        double xy=0;
        double yz=0;
        double zx=0;
        double mx=0;
        double my=0;
        double mz=0;
        for(int j=0; j<npts; j++)
        { // each point
            double x = xpt[i].at(j);
            double y = ypt[i].at(j);
            double z = zpt[i].at(j);
        //check
        if(i==cluster_check) hrec->Fill(x, y, z);

            xx += x*x/npts;
            yy += y*y/npts;
            zz += z*z/npts;
            xy += x*y/npts;
            yz += y*z/npts;
            zx += z*x/npts;
            mx += x/npts;
            my += y/npts;
            mz += z/npts;
        } //each point   

        TMatrixDSym m(3);
        m(0,0) = xx - mx*mx; 
        m(0,1) = xy - mx*my;
        m(0,2) = zx - mz*mx;
        m(1,0) = xy - mx*my;
        m(1,1) = yy - my*my;
        m(1,2) = yz - my*mz;
        m(2,0) = zx - mz*mx;
        m(2,1) = yz - my*mz;
        m(2,2) = zz - mz*mz;

        TMatrixDSymEigen me(m);
        auto& eigenval = me.GetEigenValues();
        auto& eigenvec = me.GetEigenVectors();

        double maxeval =  eigenval(0);
        int maxevalaxis = 0;
        // check
        //cout<<"axis 0 eigenvalue: "<<maxeval<<endl;

        for(int k=1; k<3; k++)
        {
            if(eigenval(k)>maxeval)
            {
                maxevalaxis = k;
                maxeval = eigenval(k);
            }
            //cout<<"axis "<<i<<" eigenvalue: "<<eigenval(i)<<endl;
        }
        
        // PCA main direction
        double pcadir_x = eigenvec(0, maxevalaxis);
        double pcadir_y = eigenvec(1, maxevalaxis);
        double pcadir_z = eigenvec(2, maxevalaxis);

        // track direction
        TVector3 pca_dir(pcadir_x, pcadir_y, pcadir_z);
        TVector3 track_dir;
        track_dir = (1./pca_dir.Mag())*pca_dir;
        
        // theta_y: 0-90 degree, upgoing
        if(track_dir.Y()<0) track_dir *= -1.0;
        
        // check
        //if( (int)(track_dir.Mag()) != 1 ) cout<<"PCA Direction vector not normalized! [units?]"<<endl;

        // Loop over point cloud to find the edges (start, end) along PCA main direction
        
        TVector3 start(xpt[i].at(0), ypt[i].at(0), zpt[i].at(0));
        TVector3 end(start);
        double start_proj = track_dir.Dot(start);
        double end_proj = track_dir.Dot(end);


        for(int j=1; j<npts; j++)
        {//each point
            TVector3 point(xpt[i].at(j), ypt[i].at(j), zpt[i].at(j));
            double point_proj = track_dir.Dot(point);
            if(point_proj < start_proj)
            {
                start=point;
                start_proj = point_proj;
            }
            if(point_proj > end_proj)
            {
                end=point;
                end_proj = point_proj;
            } 
        }//each point
        
        double costheta_y = track_dir.Y()/track_dir.Mag();
        double phi = TMath::ACos(track_dir.Z()/TMath::Sqrt(track_dir.X()*track_dir.X()+track_dir.Z()*track_dir.Z()));
        TVector3 center = 0.5*(start+end);
        double length = (end - start).Mag();
        TVector3 center2(mx, my, mz);
        TVector3 start2 = center2-0.5*length*track_dir;
        TVector3 end2 = center2+0.5*length*track_dir;

        // check
        if(i==cluster_check || cluster_check==-1){
        //cout<<"Direction: "<<track_dir.X()<<" "<<track_dir.Y()<<" "<<track_dir.Z()<<endl;
        cout<<"Cluster_id:"<<i<<endl;
        cout<<"Length: "<<length<<endl;
        cout<<"Direction: "<<costheta_y<<" "<<phi<<endl;
        cout<<"Center1: "<<center.X()<<" "<<center.Y()<<" "<<center.Z()<<endl;
        cout<<"Center2: "<<mx<<" "<<my<<" "<<mz<<endl;
        cout<<"Start1: "<<start.X()<<" "<<start.Y()<<" "<<start.Z()<<endl;
        cout<<"End1: "<<end.X()<<" "<<end.Y()<<" "<<end.Z()<<endl;
        cout<<"Start2: "<<start2.X()<<" "<<start2.Y()<<" "<<start2.Z()<<endl;
        cout<<"End2: "<<end2.X()<<" "<<end2.Y()<<" "<<end2.Z()<<endl;
        cout<<"Start Difference: "<<100*(start2.X()-start.X())/start2.X()<<" "<<100*(start2.Y()-start.Y())/start2.Y()<<" "<<100*(start2.Z()-start.Z())/start2.Z()<<endl;
        cout<<"End Difference: "<<100*(end2.X()-end.X())/end2.X()<<" "<<100*(end2.Y()-end.Y())/end2.Y()<<" "<<100*(end2.Z()-end.Z())/end2.Z()<<endl;

        leval1->SetPoint(0, start.X(), start.Y(), start.Z());
        leval1->SetPoint(1, end.X(), end.Y(), end.Z());
        leval2->SetPoint(0, start2.X(), start2.Y(), start2.Z());
        leval2->SetPoint(1, end2.X(), end2.Y(), end2.Z());
        }

    }// each cluster


    //check
    if(cluster_check!=-1){
    TCanvas* canv = new TCanvas("canv","",400,1600);
    hrec->Draw();
    leval1->Draw("same l");
    leval2->Draw("same l");
    leval1->SetLineColor(kRed);
    leval1->SetLineWidth(2);
    leval2->SetLineColor(kGreen);
    leval2->SetLineStyle(kDashed);
    leval2->SetLineWidth(2);
    canv->SaveAs("check.root");
    }

    return 0;
}
