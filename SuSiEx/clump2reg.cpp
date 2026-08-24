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
	ifstream fpi(argv[1]);
	ifstream fpref(argv[2]);
	long buff(atol(argv[3]));
	ofstream fpo(argv[4]);

	map<string, int> id2pos;
	map<string, int> id2chr;

	string line;

	while(getline(fpref, line))
	{
		istringstream cl(line);
		if(line[0] < '0' || line[0] > '9')
		{
			cerr << "Error: Please only use number to code chromosome in the bim file" << endl;
			cerr << line << endl;
			return 1;
		}
		string tid;
		int tchr, tpos;
		cl >> tchr >> tid >> line >> tpos;
		id2pos[tid] = tpos;
		id2chr[tid] = tchr;
	}

	getline(fpi, line);
	vector<bound> dat;
	while(getline(fpi, line))
	{
		vector<bound> cur;
		istringstream cl(line);
		string var, lead;
		cl >> line >> line >> lead;
		if(line == "")
			break;
		bound t;
		if(id2chr.count(lead) == 0)
		{
			cerr << "Error: Please make sure the bim file contains the lead variant: " << lead << endl;
			return 1;
		}
		t.chr = id2chr[lead];
		t.pos = id2pos[lead];
		cur.push_back(t);
		for(int j = 0 ; j < 9 ; ++j)
			cl >> line;
		if(line != "NONE")
		{
			istringstream tcl(line);
			while(getline(tcl, var, '('))
			{
				if(id2chr.count(var) == 0)
				{
					cerr << "Error: Please make sure the bim file contains the tagging variant: " << var << endl;
					return 1;
				}
				t.chr = id2chr[var];
				t.pos = id2pos[var];
				cur.push_back(t);
				getline(tcl, var, ',');
			}
		}
		sort(cur.begin(), cur.end());
		//Downstream part of the clumpped region overlapped with MHC
		if(t.chr == MHC_CHR && cur.front().pos < MHC_START && cur.back().pos > MHC_START && cur.back().pos < MHC_END)
		{
			t = cur.front();
			t.left = 1;
			t.pos -= buff;
			dat.push_back(t);
			t.left = 0;
			t.pos = MHC_START;
			dat.push_back(t);
		}
		//Upstream part of the clumpped region overlapped with MHC
		else if(t.chr == MHC_CHR && cur.front().pos < MHC_END && cur.back().pos > MHC_END && cur.front().pos > MHC_START)
		{
			t.left = 1;
			t.pos = MHC_END;
			dat.push_back(t);
			t = cur.back();
			t.left = 0;
			t.pos += buff;
			dat.push_back(t);
		}
		//The clumpped region covered the MHC
		else if(t.chr == MHC_CHR && cur.front().pos < MHC_START && cur.back().pos > MHC_END)
		{
			t = cur.front();
			t.left = 1;
			t.pos -= buff;
			dat.push_back(t);
			t.left = 0;
			t.pos = MHC_START;
			dat.push_back(t);
			t.left = 1;
			t.pos = MHC_END;
			dat.push_back(t);
			t = cur.back();
			t.left = 0;
			t.pos += buff;
			dat.push_back(t);
		}
		else
		{
			t = cur.front();
			t.left = 1;
			t.pos -= buff;
			dat.push_back(t);
			t = cur.back();
			t.left = 0;
			t.pos += buff;
			dat.push_back(t);
		}
	}
	sort(dat.begin(), dat.end());
	fpo << "Reg\tChr\tstart\tend\tstart.buffer\tend.buffer\tLength" << endl;
	int count(1);
	int n(0);
	long start(dat[0].pos);
	vector<long> len;
	for(int i = 0 ; i < dat.size(); ++i)
	{
		if(n == 0)
			start = dat[i].pos;
		if(dat[i].left)
			++n;
		else
			--n;
		if(n == 0)
		{
			fpo << count << '\t' << dat[i].chr << '\t' << start + buff << '\t' << dat[i].pos - buff << '\t' << start << '\t' << dat[i].pos << '\t' << dat[i].pos - start << endl;
			++count;
		}
	}
}
