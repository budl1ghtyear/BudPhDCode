#include "Ravl/Array1d.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/DList.hh"

using namespace RavlN;

int main(int nargs,char **argv) 
{
	Array1dC<UIntT> arr(3); arr.Fill(22);
	SArray1dC<UIntT> sarr(3); sarr.Fill(33);
	DListC<UIntT> dlist;
	dlist.Append(44);dlist.Append(44);dlist.Append(44);
	
	cout<<"Array - "<<arr<<"\n SArray - "<<sarr<<"\n DList - "<<dlist<<endl;
	
	
	
}
