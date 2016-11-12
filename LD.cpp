#include <stdio.h>
#include <string.h>     /* memcpy */
#include <stdlib.h>		/* srand, rand */
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

//GLOBALS
int numind,numlines;
int k;
char j;
char * outprefix;
string line,line1,out1,out2;
int * data;
float * pos;
float minFrac;
int * r2dist;
int bins;
//float WinSizes[10]={1000.0,5000.0,10000.0,50000.0,100000.0,500000.0,1000000.0,5000000.0,10000000.0,1000000000.0};
//float WinSizes[2]={1000.0,5000.0};
float WinSizes[10]={1000.0,2000.0,3000.0,4000.0,5000.0,10000.0,50000.0,100000.0,500000.0,1000000.0};
int main(int argc, char* argv[])
{//READ IN ARGUMENTS AND PARSE FILES
  //Check number of parameters
  if(argc != 8){
    std::cout<<"Usage: " <<argv[0]<<" AlleleCountFile" << " AllelePosFile"<<" NumberIndividuals"<<" NumberSites"<<" NumberOfBinsForR2Dist"<<" MinimumFractionGenotyped"<<" OutputPrefix"<<std::endl;
  }
  else {

    ifstream countfile (argv[1]);
    string outprefix;
    out1=string(argv[7])+".AvgR2.txt";
    out2=string(argv[7])+".Hist.txt";
    ofstream outavg (out1);
    ofstream outhist (out2);
    ifstream posfile (argv[2]);
    numind=atoi(argv[3]);
    numlines=atoi(argv[4]);
    bins=atoi(argv[5]);
    minFrac=atof(argv[6]);
    r2dist=new int[bins]();
    pos = new float[numlines*2];
    data = new int[numind*numlines]();
  if (!countfile.is_open()){
    std::cout<<"Could not open AlleleCountFile\n";
  }
  //Read data from AlleleCountFile
  else {
    int b = 0;    
    while(getline(countfile,line)){
      int a = 0;
      while(a<numind){
        data[b*numind+a]=line[2*a]-'0';
        a++;
      }      
      b++;
      }
      }
    countfile.close();
    if (!posfile.is_open()){
      std::cout<<"Could not open PosFile\n";
    }
    //Read position information from PosFile
    else {
      int c = 0;
      while(getline(posfile,line1)){
        float ps = atof(line1.c_str());
        pos[c*2]=ps;
        int zv11 = 0;
        int indcount1 = 0;
        for(int d = 0; d < numind; d++){
          if (data[c*numind+d]!=9){
            zv11=zv11+data[c*numind+d];
            indcount1++;
          }
        float m1 = float(zv11)/float(indcount1);
        pos[c*2+1]=m1;        
      }
      c++;
    }
  }
    //Cycle over window sizes specified in WinSize vector
    for (int lb=0;lb<5;lb++){
      double R2 = 0.0;
      double R22 = 0.0;
      int numobs=0;
    //Calculate r2 for given window size
      for(int p1 = 0; p1<numlines; p1++){
        for(int p2 =p1+1; p2<numlines;p2++){

  //      cout<<pos[p2*2] - pos[p1*2]<<"\n";
          if (pos[p2*2] - pos[p1*2] < WinSizes[lb+1] && pos[p2*2] - pos[p1*2] > WinSizes[lb]){
//          cout<<"here\n";
          int indcount2 = 0;
          int indcount3 = 0;
          int indcount4 = 0;
          double Zcv = 0.0;
          double Zv1 = 0.0;
          double Zv2 = 0.0;
          double r = 0.0;
          double r2 = 0.0;
          for(int ni = 0; ni<numind;ni++){
 //           cout<<data[p2*numind+ni]<<"\t"<<pos[p2*2+1]<<"\t"<<data[p1*numind+ni]<<"\t"<<pos[p1*2+1]<<"\t"<<p1<<"\t"<<p2<<"\t"<<ni<<"\n";
            if ((data[p2*numind+ni]!=9) && (data[p1*numind+ni]!=9)){
   //           cout<<"here\n";
              Zcv+=(data[p2*numind+ni]-pos[p2*2+1])*(data[p1*numind+ni]-pos[p1*2+1]);
              indcount2++;
              Zv1+=(data[p1*numind+ni]-pos[p1*2+1])*(data[p1*numind+ni]-pos[p1*2+1]);
              indcount3++;
              Zv2+=(data[p2*numind+ni]-pos[p2*2+1])*(data[p2*numind+ni]-pos[p2*2+1]);
              indcount4++;
            }
          }
  //        cout<<Zcv<<"\t"<<Zv1<<"\t"<<Zv2<<"\t"<<(float)indcount2/(float)numind<<"\t"<<minFrac<<"\t"<<indcount4<<"\n";
          
          if ((float)indcount2/(float)numind>minFrac){
 //           cout<<"here\n";
          double zcv=(double)Zcv/((double)indcount2-1.0);
          double zv1=(double)Zv1/((double)indcount3-1.0);
          double zv2=(double)Zv2/((double)indcount4-1.0);
          r=zcv/(pow(zv1,0.5)*pow(zv2,0.5));
          r2=pow(r,2);
  //        cout<<Zcv<<"\t"<<Zv1<<"\t"<<Zv2<<"\t"<<zcv<<"\t"<<zv1<<"\t"<<zv2<<"\t"<<(float)indcount2/(float)numind<<"\t"<<minFrac<<"\t"<<indcount4<<"\n";
          R2+=r2;
          R22+=pow(r2,2);
          numobs++;
          }
          float bound = 0.0;
          //Add observation to bin of r2 histogram
          for (int aa=0;aa<bins;aa++){
            if (r2>bound && r2<bound+(1.0/(float)bins)){

              r2dist[aa]+=1;
            }
            bound+=(1.0/(float)bins);
          } 
        }
      }
      if (p1%10000==0){
        cout<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<p1<<"\n";
      }
    }
    R2=R2/(double)numobs;
    R22=R22/(double)numobs;
    cout<<"Average r2 for SNPs "<<WinSizes[lb]<<" - "<<WinSizes[lb+1]<<" = "<<R2<<"\n";
    cout<<"Number of comparisons = "<<numobs<<"\n";
    outavg<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<R2<<"\t"<<R22-pow(R2,2)<<"\n";
    float bound1 = 0.0;
    for (int cc=0;cc<bins;cc++){
      outhist<<WinSizes[lb]<<"\t"<<WinSizes[lb+1]<<"\t"<<bound1<<"\t"<<bound1+(1.0/(float)bins)<<"\t"<<r2dist[cc]<<"\n";
      bound1+=(1.0/(float)bins);
    }
  }
  }
  return 0;
}

//    for (int g = 0; g<numlines*2;g++){
//      cout<<pos[g]<<"\n";
//    }
    
//    cout<<"WindowSize = "<<WSLB<<" - "<<WSUB<<"\n";
//    cout<<"Covariance = "<<cv <<"\n";
//    cout<<"V1 = "<<v1<<"\n";
//    cout<<"V2 = "<<v2<<"\n";
//    cout<<"r = "<<r<<"\n";
//    cout<<"r2 = "<<r2<<"\n";  

