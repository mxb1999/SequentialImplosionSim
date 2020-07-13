#include <fstream>
#include <iostream>

using namespace std;

int determineType(string type);
bool checkVals(string in, void* arr, int type, int* loclist);
void checkArr(void* arr, int* dimlist, int dimsize, string type, string file, string name);
