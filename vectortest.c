int matrixMult(int a[], int b[])
{
  int N = sizeof(a)/sizeof(a[0]);
for(int i=0;i<N;i+=2)
{
  a[i::i+2] = a[i::i+2] * b[i::i+2]
}
}
int main(int argc,char** argv)
{
  int a[] = {1,3,5,6};
  int b[] = {1,3,5,6};
  printf("%s%i",matrixMult(a,b));
}
