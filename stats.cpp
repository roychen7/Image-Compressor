
#include "stats.h"

stats::stats(PNG & im){
  int width = im.width();
  int height = im.height();
  sumHueX = vector<vector<double>>(height,vector<double>(width,0));
  sumHueY = vector<vector<double>>(height,vector<double>(width,0));
  sumSat = vector<vector<double>>(height,vector<double>(width,0));
  sumLum = vector<vector<double>>(height,vector<double>(width,0));
  hist = vector<vector<vector<int>>>(height,vector<vector<int>>(width,vector<int>(36,0)));
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      sumHueX[i][j] = cos(im.getPixel(j,i)->h * PI/180);
      sumHueY[i][j] = sin(im.getPixel(j,i)->h * PI/180);
      sumSat[i][j] = im.getPixel(j,i)->s;
      sumLum[i][j] = im.getPixel(j,i)->l;
      if((i == 0) && (j != 0)){
        sumHueX[i][j] += sumHueX[i][j-1];
        sumHueY[i][j] += sumHueY[i][j-1];
        sumSat[i][j] += sumSat[i][j-1];
        sumLum[i][j] += sumLum[i][j-1];
        hist[i][j] = hist[i][j-1];
      } else if((j == 0) && (i != 0)){
        sumHueX[i][j] += sumHueX[i-1][j];
        sumHueY[i][j] += sumHueY[i-1][j];
        sumSat[i][j] += sumSat[i-1][j];
        sumLum[i][j] += sumLum[i-1][j];
        hist[i][j] = hist[i-1][j];
      } else if((j != 0) && (i != 0)){
        sumHueX[i][j] += sumHueX[i-1][j];
        sumHueX[i][j] += sumHueX[i][j-1];
        sumHueX[i][j] -= sumHueX[i-1][j-1];
        sumHueY[i][j] += sumHueY[i-1][j];
        sumHueY[i][j] += sumHueY[i][j-1];
        sumHueY[i][j] -= sumHueY[i-1][j-1];
        sumSat[i][j] += sumSat[i-1][j];
        sumSat[i][j] += sumSat[i][j-1];
        sumSat[i][j] -= sumSat[i-1][j-1];
        sumLum[i][j] += sumLum[i-1][j];
        sumLum[i][j] += sumLum[i][j-1];
        sumLum[i][j] -= sumLum[i-1][j-1];
        for(int k = 0; k < 36; k++){
          hist[i][j][k] += hist[i-1][j][k];
          hist[i][j][k] += hist[i][j-1][k];
          hist[i][j][k] -= hist[i-1][j-1][k];
        }
      }
      for(int k = 0; k < 36; k++){
        if((k*10 <= im.getPixel(j,i)->h) && (im.getPixel(j,i)->h < (k+1)*10)){
          hist[i][j][k]++;
        }
      }
    }
  }
}

long stats::rectArea(pair<int,int> ul, pair<int,int> lr){
  int x1 = ul.first;
  int x2 = lr.first;
  int y1 = ul.second;
  int y2 = lr.second;
  int width = x2-x1+1;
  int length = y2-y1+1;
  return abs(width*length);

}

HSLAPixel stats::getAvg(pair<int,int> ul, pair<int,int> lr){
  int x1 = ul.first;
  int x2 = lr.first;
  int y1 = ul.second;
  int y2 = lr.second;
  double avgSat, avgLum, avgHueX, avgHueY;
  // Four cases: Top, Left, Corner, arbitrary.
  long numPix = rectArea(ul,lr);
  if((x1 == 0) && (y1 == 0)){
    avgSat = sumSat[y2][x2]/numPix;
    avgLum = sumLum[y2][x2]/numPix;
    avgHueX = sumHueX[y2][x2]/numPix;
    avgHueY = sumHueY[y2][x2]/numPix;
  } else if(x1 == 0){
    avgSat = (sumSat[y2][x2]-sumSat[y1-1][x2])/numPix;
    avgLum = (sumLum[y2][x2]-sumLum[y1-1][x2])/numPix;
    avgHueX = (sumHueX[y2][x2]-sumHueX[y1-1][x2])/numPix;
    avgHueY = (sumHueY[y2][x2]-sumHueY[y1-1][x2])/numPix;
  } else if(y1 == 0){
    avgSat = (sumSat[y2][x2]-sumSat[y2][x1-1])/numPix;
    avgLum = (sumLum[y2][x2]-sumLum[y2][x1-1])/numPix;
    avgHueX = (sumHueX[y2][x2]-sumHueX[y2][x1-1])/numPix;
    avgHueY = (sumHueY[y2][x2]-sumHueY[y2][x1-1])/numPix;
  } else{
  avgSat = (sumSat[y2][x2] - sumSat[y2][x1-1]
                  - sumSat[y1-1][x2] + sumSat[y1-1][x1-1])/numPix;
  avgLum = (sumLum[y2][x2] - sumLum[y2][x1-1]
                  - sumLum[y1-1][x2] + sumLum[y1-1][x1-1])/numPix;
  avgHueX = (sumHueX[y2][x2] - sumHueX[y2][x1-1]
                  - sumHueX[y1-1][x2] + sumHueX[y1-1][x1-1])/numPix;
  avgHueY = (sumHueY[y2][x2] - sumHueY[y2][x1-1]
                  - sumHueY[y1-1][x2] + sumHueY[y1-1][x1-1])/numPix;
                }
  double avgHue = atan2(avgHueY,avgHueX) * 180/PI;
  if (avgHue < 0){
    avgHue += 360;
  }
  return HSLAPixel(avgHue,avgSat,avgLum,1.0);
}

vector<int> stats::buildHist(pair<int,int> ul, pair<int,int> lr){
  int x1 = ul.first;
  int x2 = lr.first;
  int y1 = ul.second;
  int y2 = lr.second;
  vector<int> output(36,0);
  if((x1 == 0) && (y1 == 0)){
    for(int k = 0; k < 36; k++){
      output[k] += hist[y2][x2][k];
    }
  } else if(x1 == 0){
    for(int k = 0; k < 36; k++){
      output[k] += hist[y2][x2][k];
      output[k] -= hist[y1-1][x2][k];
    }
  } else if(y1 == 0){
    for(int k = 0; k < 36; k++){
      output[k] += hist[y2][x2][k];
      output[k] -= hist[y2][x1-1][k];
    }
  } else{
    for(int k = 0; k < 36; k++){
      output[k] += hist[y2][x2][k];
      output[k] -= hist[y2][x1-1][k];
      output[k] -= hist[y1-1][x2][k];
      output[k] += hist[y1-1][x1-1][k];
    }
  }
  return output;
}

// takes a distribution and returns entropy
// partially implemented so as to avoid rounding issues.
double stats::entropy(vector<int> & distn,int area){
    double entropy = 0.;
    for (int i = 0; i < 36; i++) {
        if (distn[i] > 0 )
            entropy += ((double) distn[i]/(double) area)
                                    * log2((double) distn[i]/(double) area);
    }
    return  -1 * entropy;

}

double stats::entropy(pair<int,int> ul, pair<int,int> lr){
  int area = rectArea(ul,lr);
  vector<int> entropyHist = buildHist(ul,lr);
  return entropy(entropyHist, area);
}
