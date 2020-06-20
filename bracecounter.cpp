#include <iostream>
#include <fstream>
using namespace std;
int main(int argc, char const *argv[]) {
  ifstream txt ("implSim.txt");
  int fcount = 0;
  int bcount = 0;
  while(txt)
  {
    string input;
    txt >> input;
    int s = input.size();
    for(int i = 0; i < s; i++)
    {
      char curr = input.at(i);
      if(curr == '{')
      {
        fcount++;
      }else if(curr == '}')
      {
        bcount++;
      }
    }
  }
  cout << "{ Count: " << fcount<<endl;
  cout << "} Count: " << bcount<<endl;
  return 0;
}
