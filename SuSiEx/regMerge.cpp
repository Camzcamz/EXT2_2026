# include <iostream>
# include <fstream>
# include <sstream>
# include <vector>
# include <algorithm>
# include <map>
# include <bitset>

using namespace std;

const int MHC_CHR = 6;
const long MHC_START = 28477797;
const long MHC_END = 33448354;

class bound
{
	public:
		int chr;
		long pos;
		int left;
		int dat;
};

bool operator<(const bound &a, const bound& b)
{
	if(a.chr != b.chr)
		return a.chr < b.chr;
	if(a.pos != b.pos)
		return a.pos < b.pos;
	return a.left > b.left;
}

int main(int argc, char **argv)
{
	ifstream fplist(argv[1]);
	ofstream fpo(argv[2]);
	
	vector<string> path, lab;
	string tpath, tlab;
	while(fplist >> tpath >> tlab)
	{
		path.push_back(tpath);
		lab.push_back(tlab);
	}
	
	vector<bound> dat;

	for(int i = 0 ; i < path.size(); ++i)
	{
		ifstream fpi(path[i].c_str());
		if(!fpi)
		{
			cerr << "Error: cannot open " << path[i] << endl;
			return 1;
		}
		string line;
		getline(fpi, line);
		bound a, b;
		a.left = 1;
		b.left = 0;
		a.dat = i;
		b.dat = i;
		while(getline(fpi, line))
		{
			istringstream cl(line);
			cl >> line >> a.chr >> line >> line >> a.pos >> b.pos;
			b.chr = a.chr;
			dat.push_back(a);
			dat.push_back(b);
		}
	}
	sort(dat.begin(), dat.end());
	fpo << "Reg\tChr\tstart.buffer\tend.buffer";
	for(int i = 0 ; i < lab.size(); ++i)
		fpo << '\t' << lab[i] << ".length";
	fpo << endl;
	int count(1);
	long start(dat[0].pos);
	vector<long> len;
	vector<int> n;
	len.resize(lab.size(), 0);
	n.resize(lab.size(), 0);
	int alln(0);
	for(int i = 0 ; i < dat.size(); ++i)
	{
		if(alln == 0)
			start = dat[i].pos;
		for(int j = 0 ; j < lab.size(); ++j)
		{
			if(n[j])
				len[j] += dat[i].pos - dat[i - 1].pos;
		}
		if(dat[i].left)
		{
			++n[dat[i].dat];
			++alln;
		}
		else
		{
			--n[dat[i].dat];
			--alln;
		}
		
		if(alln == 0)
		{
			fpo << count << '\t' << dat[i].chr << '\t' << start << '\t' << dat[i].pos;
			for(int j = 0 ; j < lab.size(); ++j)
				fpo << '\t' << len[j];
			fpo << endl;
			for(int j = 0 ; j < lab.size(); ++j)
				len[j] = 0;
			++count;
		}
	}
}
