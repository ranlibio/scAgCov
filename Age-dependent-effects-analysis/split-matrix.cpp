#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <math.h>

using namespace std;

void parseStr(string& ss, char c, vector<string>& vc2)
{
	vc2.clear();
	for (int i = 0;; i++)
	{
		string::size_type pos = ss.find(c);
		if (pos == string::npos)
		{
			vc2.push_back(ss);
			break;
		}
		else
		{
			string s = ss.substr(0, pos);
			string s1 = ss.substr(pos + 1, ss.size() - pos - 1);
			vc2.push_back(s);
			ss = s1;
		}
	}
}

void split(string& filename, vector<string>& barcode)
{
	string task = filename;
	ifstream fin(task.c_str());
	vector<string> vs;
	while (!fin.eof())
	{
		char s[1000];
		fin.getline(s, 1000);
		string ss = s;
		if (ss.size() == 0) break;
		parseStr(ss, ',', vs);
		barcode.push_back(vs[0]);
	}
}

int main(int argc, char* argv[])
{
	string task1 = argv[1];
	string task2 = argv[2];
	string task3 = argv[3];
	string task4 = argv[4];

	vector<string> barcode;
	split(task1,barcode);
	
	int start=atoi(barcode[0].c_str());
	int end=atoi(barcode[barcode.size()-1].c_str());
	int len=end-start+1;
	
	cout<<start<<"\t"<<end<<"\t"<<len<<endl;

	ifstream fin1(task2.c_str());
	ofstream fout1(task3.c_str());

	vector<string> vs;
	double line = 0;
	int cellCount=0;
	double step=100000;
	double skip=4;
	while (!fin1.eof())
	{
		++line;
		int remainder = fmod(line,100000);
		if (remainder == 0) cout << line << endl;
		char s[200];
		fin1.getline(s, 200);
		string ss = s;
		if(line<skip) continue;
		if (ss.size() == 0) break;
		parseStr(ss ,' ', vs);
		
		int location=atoi(vs[1].c_str());
		if(location>end) break;
		    
		if(location>=start && location<=end)
		{
			++cellCount;
			fout1 << vs[0] << "\t" << location-start+1 << "\t" << vs[2] << endl;
		}
	}
	fin1.close();
	fout1.close();

	
	ifstream fin2(task3.c_str());
	ofstream fout2(task4.c_str());
	int k=0;
	while (!fin2.eof())
	{
		++k;
		char s[200];
		fin2.getline(s, 200);
		string ss = s;
		if (ss.size() == 0) break;
		
		if(k==1)
		{
			fout2<<"%%MatrixMarket matrix coordinate real general"<<"\n"<<"%"<<endl;
			fout2 << 27943 << "\t" << len <<"\t"<< cellCount << endl;
			fout2<<ss<<endl;
		}else{
			fout2<<ss<<endl;
		}
	}
	fin2.close();
	fout2.close();
}

