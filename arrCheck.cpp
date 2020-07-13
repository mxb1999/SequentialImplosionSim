#include <iostream>
#include <fstream>
#include "arrCheck.h"
#include <math.h>

int determineType(string type)
{
  if(type.compare("double*") == 0)
  {
    return 0;
  }else if(type.compare("double**") == 0)
  {
    return 1;
  }else if(type.compare("double***") == 0)
  {
    return 2;
  }else if(type.compare("double****") == 0)
  {
    return 3;
  }else if(type.compare("int*") == 0)
  {
    return 4;
  }else if(type.compare("int**") == 0)
  {
    return 5;
  }else if(type.compare("int***") == 0)
  {
    return 6;
  }else if(type.compare("int****") == 0)
  {
    return 7;
  }else
  {
    return -1;
  }
};
bool checkVals(string in, void* arr, int type, int* loclist)
{
  double epsilon = 1e-4;
  switch(type)
  {
    case 0:{
      double* arrDub = (double*)arr;
      double toCheck = arrDub[loclist[0]];
      if(fabs(stod(in)-toCheck) <= epsilon * fabs(stod(in)))
      {
        return true;
      }else
      {
        printf("%s%f%s%f\n", "Failed: ", toCheck, " || Actual: ", stod(in));
        return false;
      }}
    case 1:{
      double** arrDub = (double**)arr;
      double toCheck = arrDub[loclist[0]][loclist[1]];
      if(fabs(stod(in)-toCheck) <= epsilon * fabs(stod(in)))
      {
        return true;
      }else
      {
        printf("%s%f%s%f\n", "Failed: ", toCheck, " || Actual: ", stod(in));
        return false;
      }}
    case 2:
      {
      double*** arrDub = (double***)arr;
      double toCheck = arrDub[loclist[0]][loclist[1]][loclist[2]];
      if(fabs(stod(in)-toCheck) <= epsilon * fabs(stod(in)))
      {
        return true;
      }else
      {
        printf("%s%f%s%f\n", "Failed: ", toCheck, " || Actual: ", stod(in));
        return false;
      }}
    case 3:{
      double**** arrDub = (double****)arr;
      double toCheck = arrDub[loclist[0]][loclist[1]][loclist[2]][loclist[3]];
      if(fabs(stod(in)-toCheck) <= epsilon * fabs(stod(in)))
      {

        return true;
      }else
      {
        printf("%s%f%s%f\n", "Failed: ", toCheck, " || Actual: ", stod(in));
        return false;
      }}
    case 4:
      {
      int* intarr = (int*)arr;
      int toCheck = intarr[loclist[0]];
      if(stoi(in) == toCheck)
      {
        return true;
      }else
      {
        printf("%s%d%s%d\n", "Failed: ", toCheck, " || Actual: ", stoi(in));
        return false;
      }}
    case 5:
      {
      int** intarr = (int**)arr;
      int toCheck = intarr[loclist[0]][loclist[1]];
      if(stoi(in) == toCheck)
      {
        return true;
      }else
      {
        printf("%s%d%s%d\n", "Failed: ", toCheck, " || Actual: ", stoi(in));
        return false;
      }}
    case 6:
      {
      int*** intarr = (int***)arr;
      int toCheck = intarr[loclist[0]][loclist[1]][loclist[2]];
      if(stoi(in) == toCheck)
      {
        return true;
      }else
      {
        printf("%s%d%s%d\n", "Failed: ", toCheck, " || Actual: ", stoi(in));
        return false;
      }}
    case 7:
      {
      int**** intarr = (int****)arr;
      int toCheck = intarr[loclist[0]][loclist[1]][loclist[2]][loclist[3]];
      if(stoi(in) == toCheck)
      {
        return true;
      }else
      {
        printf("%s%d%s%d\n", "Failed: ", toCheck, " || Actual: ", stoi(in));
        return false;
      }
    }
  }
  return false;
};
void checkArr(void* arr, int* dimlist, int dimsize, string type, string file, string name)
{
  bool hasReached = false;
  ifstream input (file);
  ofstream output ("output.txt");
  int arrType = determineType(type);
  string in;
  int track = 0;
  while(getline(input, in))
  {
    if(in.find_last_of(name) != string::npos)
    {

      hasReached = true;
      int loc = in.find_last_of(" ");
      int locarr[dimsize] = {0};
      int trackclone = track;
      int dimprod = 1;
      for(int i = 0; i < dimsize - 1; i++)
      {
        dimprod *= dimlist[i];
      }
      for(int i = 0; i < dimsize; i++)
      {
        if(i == dimsize - 1)
        {
          locarr[i] = trackclone % dimlist[i];
        }else
        {
          locarr[i] = trackclone/dimprod;
          dimprod /= dimlist[i+1];
        }
      }
      if(!checkVals(in.substr(loc+1), arr, arrType, locarr))
      {
        output << track <<"\n";
      }
    } else if(hasReached)
    {
      break;
    }
    track++;
  }
  input.close();
  output.close();

};
